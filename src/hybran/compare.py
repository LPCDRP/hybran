import itertools
import logging
import os
import sys

from intervaltree import Interval, IntervalTree
import networkx as nx

from . import (
    designator,
    extractor,
    fileManager,
)
from .bio import SeqIO


def main(args):

    logger = logging.getLogger('Compare')

    gbk_file1, gbk_file2 = args.annotations
    outdir = args.outdir

    # We preserve the original basename here in case in the future we decide to support matching contigs for comparison by their names.
    # In that case, the file name prefix might be part of the sequence ID that we can use to help with match-making.
    og_basename1, ext1 = os.path.splitext(os.path.basename(gbk_file1))
    og_basename2, ext2 = os.path.splitext(os.path.basename(gbk_file2))

    if (
            ext1[1:] not in fileManager.exts['genbank'] or
            ext2[1:] not in fileManager.exts['genbank']
    ):
        sys.exit(f"Input file must be in Genbank format")
    if og_basename1 == og_basename2:
        prefix1 = os.path.splitext(gbk_file1)[0]
        prefix2 = os.path.splitext(gbk_file2)[0]
        # same file (or only differing in extension. like x.gbk vs x.gb)
        if prefix1 == prefix2:
            basename1 += "_1"
            basename2 += "_2"
        # same file name but in different folders.
        # This is a common scenario for different version runs of the same sample.
        else:
            basename1 = prefix1.replace(os.sep, '_')
            basename2 = prefix2.replace(os.sep, '_')
    else:
        basename1 = og_basename1
        basename2 = og_basename2

    record_list = generate_record(gbk_file1)
    alt_record_list = generate_record(gbk_file2)

    n_records = len(record_list)
    n_alt_records = len(alt_record_list)
    if n_records != n_alt_records:
        min_n = min([n_records, n_alt_records])
        logger.warning((
            f"{basename1} has {n_records} sequences and {basename2} has {n_alt_records}. "
            f"Only comparing the first {min_n} to each other..."
        ))
        n_records = min_n

    matching = []
    conflicts = []
    uniques = []
    alt_uniques = []
    ovl_graph = nx.Graph()
    for chrom in range(n_records):
        (
            chrom_matching,
            chrom_conflicts,
            chrom_nonconfl,
            chrom_alt_nonconfl,
            chrom_uniques,
            chrom_alt_uniques,
            chrom_ovl_graph,
        ) = compare(
            record_list[chrom],
            alt_record_list[chrom],
        )
        matching.extend(chrom_matching)
        conflicts.extend(chrom_conflicts)
        uniques += chrom_nonconfl + chrom_uniques
        alt_uniques += chrom_alt_nonconfl + chrom_alt_uniques
        ovl_graph = nx.compose(ovl_graph, chrom_ovl_graph)

    feature_list = list(itertools.chain.from_iterable(record_list))
    alt_feature_list = list(itertools.chain.from_iterable(alt_record_list))

    write_reports(
        matching,
        conflicts,
        uniques,
        alt_uniques,
        ovl_graph,
        feature_list,
        alt_feature_list,
        basename1,
        basename2,
        outdir,
    )

def complementary(f1, f2, f1_status=None, f2_status=None):
    """
    Given two (Autarkic)SeqFeatures, check whether their reference-correspondence is complementary.
    That is, one has a good start and bad stop while the other has a bad start and good stop.
    :param f1: SeqFeature
    :param f2: SeqFeature
    :param f1_status: bool 2-tuple (goodstart, goodstop) for f1
    :param f2_status: bool 2-tuple (goodstart, goodstop) for f2
        If None, attempt to use (f1(2).rcs, f1(2).rce). The argument allows flexibility for
        providing coord_check results for f1 or f2 with respect to alternative reference genes.
    :return: bool whether the two features are complementary
    """
    if f1_status is None:
        f1_status = (f1.rcs, f1.rce)
    if f2_status is None:
        f2_status = (f2.rcs, f2.rce)
    return (
        (any(f1_status) and any(f2_status))
        and (int(f1_status[0])+int(f2_status[0]), int(f1_status[1])+int(f2_status[1]))==(1,1)
    )

def overlap_inframe(loc1, loc2):
    """
    Say whether two FeatureLocations are overlapping and share the same reading frame.
    This is for the purpose of CDSs to determine whether they should be corresponding to the
    same genomic feature or whether one gene has lost its stop codon and merged into its neighbor.

    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features overlap in-frame
    """
    # overlap determination adapted from https://stackoverflow.com/a/2953979
    def overlap_parts(loc1, loc2):
        overlap = (min(loc1.end, loc2.end) - max(loc1.start, loc2.start))

        if loc1.strand == loc2.strand and overlap > 0:
            # pseudogenes may occupy multiple reading frames, so
            # check both the start-defined and stop-defined frames,
            # as well as the actual overlapping part.
            if (loc1.start == loc2.start
                or loc1.end == loc2.end
                or (
                    ((loc1.start - loc2.start) % 3 == 0
                     or (loc1.end - loc2.end) % 3 == 0)
                    and (overlap % 3 == 0))
            ):
                return True
        return False

    for i in loc1.parts:
        for j in loc2.parts:
            if overlap_parts(i,j):
                return True
    return False

def have_same_stop(loc1, loc2):
    """
    Say whether two FeatureLocations have the same stop position.
    This is a special case of overlap_inframe().
    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features have the same stop position
    """
    if loc1.strand == loc2.strand:
        return (
            (loc1.strand == 1 and loc1.end == loc2.end)
            or
            (loc1.strand == -1 and loc1.start == loc2.start)
        )
    return False

def hybran_np(feature):
    """
    Note parser to identify types of pseudogenes specifically from Hybran annotation files
    :param feature: A SeqFeature
    :return evidence_code: String of the evidence code generated from pseudoscan
    """
    pseudo_notes = [_ for _ in feature.qualifiers['note'] if "Hybran/Pseudoscan:evidence:" in _]
    # older hybran releases don't have these evidence code notes
    if not pseudo_notes:
        return '.'
    evidence_code = pseudo_notes[0].split(":", 2)[2]
    # A long string gets line-wrapped in the genbank file and Biopython doesn't correctly undo the line-wrapping.
    evidence_code = evidence_code.replace(' ', '')
    return evidence_code

pgap_evidence_codes = [
    'too many ambiguous residues',
    'frameshifted',
    'incomplete',
    'internal stop',
]

def pgap_np(feature):
    """
    Note parser to identify types of pseudogenes specifically from PGAP annotation files
    :param feature: A SeqFeature
    :return evidence_code: String from the note field of a CDS with a '/pseudo' qualifier
    """
    evid = {}
    evidence_code = []
    for note in feature.qualifiers['note']:
        components = note.split(';')
        for component in components:
            component = component.strip()
            if component in pgap_evidence_codes:
                evidence_code.append(component.replace(' ', '_'))
            else:
                break
    evidence_code = f"{';'.join(evidence_code) if evidence_code else '.'}"
    return evidence_code

def which_np(gbk_record):
    for header, comments in gbk_record.annotations.items():
        header = str(header).lower()
        comments = str(comments).lower()
        if ('hybran' in header) or ('hybran' in comments):
            return hybran_np
        if ('pgap' in header) or ('pgap' in comments):
            return pgap_np

    #unknown or generic note parser return null entry.
    return lambda _:'.'

def generate_record(gbk):
    """
    Creates a list of of lists of features for each contig represented in a .gbk annotation file.
    :param gbk: String of a path to a .gbk annotation file.
    :return anno_list: List of lists of features from the input.gbk file ordered by position.
    """
    all_records = SeqIO.parse(gbk, "genbank")
    anno_list = []
    int_tree = IntervalTree()

    for anno_record in all_records:
        np = which_np(anno_record)
        feature_list = []
        for f in anno_record.features:
            f_interval = get_feature_interval(f)
            if f.type == 'CDS' and f_interval not in int_tree:
                if designator.is_pseudo(f.qualifiers):
                    f.evidence = np(f)
                feature_list.append(f)
                int_tree.add(f_interval)
        anno_list.append(feature_list)
    return anno_list

def get_pseudo_type(feature):
    return f"{feature.evidence if hasattr(feature, 'evidence') else '.'}"

def get_feature_interval(feature):
    """
    :param feature: A SeqFeature
    :return Interval object of feature
    """
    return Interval(int(feature.location.start), int(feature.location.end))

def compare(feature_list, alt_feature_list, eliminate_colocated=True):
    """
    Organize each feature into various categories based on overlaps between the annotation files.

    :param feature_list: List of features from the input.gbk file ordered by position.
    :param alt_feature_list: List of features from the alternative input.gbk file ordered by position.
    :param eliminate_colocated: bool whether to explore further conflicts with colocated genes.
      If True, you will consider colocated features as accounted for and instead see what would have been the additional conflicting features as unique or nonconflicting overlapping.
    :returns:
      - co_located: List of 2-tuples of SeqFeatures that have exactly matching positions.
      - conflicting: List of 2-tuples of SeqFeatures that conflict (for CDSs, meaning they overlap in-frame).
      - nonconfl_overlap: List of SeqFeatures from the first annotation that are putatively non-conflicting with the second
          (for CDSs, meaning they overlap out-of-frame with at least one other feature)
      - alt_nonconfl_overlap: List of SeqFeatures from the second annotation that are putatively non-conflicting with the first
      - unique: List of SeqFeatures that are unique to the first annotation file, not overlapping anything in the second.
      - alt_unique: List of SeqFeatures that are unique to the second annotation file, not overlapping anything in the first.
      - G_all: networkx unidirected bipartite graph representing all overlap relationships (excluding identical overlaps if eliminate_colocated).
    """

    # Use a bipartite graph to relate conflicts
    # We use a directed graph here because we want to be able to rely on edges being F1 -> F2
    # so that we can identify the source genome easily.
    #
    # for node attributes:
    # bipartite=0 => genome 1
    # bipartite=1 => genome 2
    G_conflict = nx.DiGraph()

    alt_interval_tree = IntervalTree()
    for i, alt_feature in enumerate(alt_feature_list):
        alt_feature.label = f"Y{i}"
        G_conflict.add_node(
            alt_feature.label,
            annotation=alt_feature,
            bipartite=1,
        )
        alt_interval_tree.addi(int(alt_feature.location.start), int(alt_feature.location.end), alt_feature)

    # A second bipartite graph to track all overlaps (except exact matches, unless not eliminate_colocated)
    # we don't want edge direction here because we want to easily identify edges originating from
    # features from either genome.
    G_all = G_conflict.copy().to_undirected()

    co_located = []
    conflicting = []
    nonconfl_overlap = []
    alt_nonconfl_overlap = []
    unique = []
    alt_unique = []

    for i, feature in enumerate(feature_list):
        feature.label = f"X{i}"
        for G in [G_conflict, G_all]:
            G.add_node(
                feature.label,
                annotation=feature,
                bipartite=0,
            )
        alt_overlapping_intervals = alt_interval_tree.overlap(get_feature_interval(feature))
        if alt_overlapping_intervals:
            alt_overlappers = sorted(
                [alt_feature for alt_start, alt_end, alt_feature in alt_overlapping_intervals],
                key=lambda f: f.location.start
            )

        colo = []
        conf = []
        non_conf = []
        for alt_f in alt_overlappers:
            if feature.location == alt_f.location:
                colo.append(alt_f)
            # we can currently only identify conflicts for CDSs
            elif (
                    feature.type == 'CDS'
                    and feature.type == alt_f.type
                    and overlap_inframe(feature.location, alt_f.location)
            ):
                conf.append(alt_f)
            else:
                non_conf.append(alt_f)

        # 'conflicting' : overlaps in frame
        for alt_feature in conf:
            for G in [G_conflict, G_all]:
                # alt_feature could have been dropped from the graph in a previous match.
                # if so, don't add it back into the graph.
                if alt_feature.label in G.nodes:
                    G.add_edge(feature.label, alt_feature.label, inframe=True, colo=False)

        # 'non-conflicting' : overlaps out of frame
        for alt_feature in non_conf:
            G_all.add_edge(feature.label, alt_feature.label, inframe=False, colo=False)

        # 'co_located' : exact match
        # looping, but the only time len(colo)>1 is if there are redundant annotation entries.
        for alt_feature in colo:
            co_located.append( (feature, alt_feature) )
            if eliminate_colocated:
                # we don't care about further overlaps with these pairs, as they're fully accounted for.
                G_conflict.remove_nodes_from([alt_feature.label, feature.label])
            else:
                for G in [G_conflict, G_all]:
                    G.add_edge(feature.label, alt_feature.label, inframe=True, colo=True)

    #
    # See what's left now that the dust has settled
    #
    for comp in nx.weakly_connected_components(G_conflict):
        comp = list(comp)
        # "connected components" with only one node represent unique features.
        if len(comp) == 1:
            annotation = G_conflict.nodes[comp[0]]['annotation']
            # If it's also zero-degree in the nonconfl overlap graph, it's truly unique: not overlapping anything
            if G_all.degree[comp[0]] == 0:
                if G_conflict.nodes[comp[0]]['bipartite'] == 0:
                    unique.append(annotation)
                else:
                    alt_unique.append(annotation)
            # putatively nonconflicting, but overlapping something
            else:
                if G_conflict.nodes[comp[0]]['bipartite'] == 0:
                    nonconfl_overlap.append(annotation)
                else:
                    alt_nonconfl_overlap.append(annotation)
        else:
            conflicting += [
                (G_conflict.nodes[x]['annotation'], G_conflict.nodes[y]['annotation'])
                for x,y in G_conflict.edges(comp) if not G_conflict[x][y]['colo']
            ]

    return (
        co_located,
        conflicting,
        nonconfl_overlap,
        alt_nonconfl_overlap,
        unique,
        alt_unique,
        G_all,
    )

#
# Helper functions for write_reports()
#

def unpack_feature(feature):
    return {
        'ltag': extractor.get_ltag(feature),
        'gene': extractor.get_gene(feature, tryhard=False),
        'start': str(feature.location.start + 1),
        'end': str(feature.location.end),
        'strand': str(feature.location.strand),
        'pseudo': str(int(designator.is_pseudo(feature.qualifiers))),
        'pseudo_type': get_pseudo_type(feature),
    }

def format_match(pair, G=None):
    f1, f2 = [unpack_feature(_) for _ in pair]
    return [
        f1['ltag'], f1['gene'], f1['pseudo'], f1['pseudo_type'],
        f2['ltag'], f2['gene'], f2['pseudo'], f2['pseudo_type'],
        f1['start'], f1['end'], f1['strand'],
        assign_tags(pair),
    ]

def format_conflict(pair, G=None):
    f1, f2 = [unpack_feature(_) for _ in pair]
    return [
        f1['ltag'], f1['gene'], f1['start'], f1['end'], f1['strand'], f1['pseudo'], f1['pseudo_type'],
        f2['ltag'], f2['gene'], f2['start'], f2['end'], f2['strand'], f2['pseudo'], f2['pseudo_type'],
        assign_tags(pair),
    ]

def assign_tags(pair):
    """
    Create a string to label special properties of a pair of opposing features to facilitate filtering.

    :param pair: a 2-tuple of SeqFeatures
    :return: str ;-delimited properties applying to a pair.
    """
    tags = []
    if name(pair[0]) and not name(pair[1]):
        tags.append("exclusively_named_by_1")
    elif not name(pair[0]) and name(pair[1]):
        tags.append("exclusively_named_by_2")
    elif name(pair[0]) is not None and name(pair[0]) == name(pair[1]):
        tags.append("identically_named")

    f1_pseudo = designator.is_pseudo(pair[0].qualifiers)
    f2_pseudo = designator.is_pseudo(pair[1].qualifiers)
    if f1_pseudo or f2_pseudo:
        tags.append("involving_pseudo")
    elif f1_pseudo and not f2_pseudo:
        tags.append("only_pseudo_in_1")
    elif not f1_pseudo and f2_pseudo:
        tags.append("only_pseudo_in_2")

    if not tags:
        return "."
    else:
        return ';'.join(tags)

def name(feature):
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    else:
        return None

def get_named(feature_list, key=lambda _:_, unique=False):
    """
    Get all features that have a name defined.

    :param feature_list: list of SeqFeatures or list of containers including SeqFeatures
    :param key: function that returns the SeqFeature from an item in feature_list
    :param unique: bool whether to return a list of unique SeqFeatures meeting the criteria
    """
    named_features = [key(_) for _ in feature_list if name(key(_))]
    if unique:
        unique_features = []
        for _ in named_features:
            if _ not in unique_features:
                unique_features.append(_)
        named_features = unique_features
    return named_features

def differentially_named(feature_tuple_list):
    """
    :param feature_tuple_list: list of 2-tuples of SeqFeatures
    :return: two lists of feature pairs.
    One list where only the first element in the tuple is named and another
    list where only the second element in the tuple is named.
    """
    f1_exclusive = []
    f2_exclusive = []
    for f1, f2 in feature_tuple_list:
        if name(f1) and not name(f2):
            f1_exclusive.append((f1, f2))
        elif not name(f1) and name(f2):
            f2_exclusive.append((f1, f2))
    return f1_exclusive, f2_exclusive

def identically_named(feature_tuple_list):
    same_names = []
    for f1, f2 in feature_tuple_list:
        if name(f1) is not None and name(f1) == name(f2):
            same_names.append((f1, f2))
    return same_names

def format_unique(feature, G):
    f = unpack_feature(feature)
    overlaps = [G.nodes[_[1]]['annotation'] for _ in G.edges(feature.label)]
    overlap_note = ";".join([
        f"{extractor.get_gene(_)}|{str(_.location)}|pseudo={int(designator.is_pseudo(_.qualifiers))}|"
        f"pseudo_type={get_pseudo_type(_)}" for _ in overlaps
    ])
    if not overlap_note:
        overlap_note = "."
    return [
        f['ltag'],
        f['gene'],
        f['start'],
        f['end'],
        f['strand'],
        f['pseudo'],
        f['pseudo_type'],
        overlap_note,
    ]


def write_reports(
        matching,
        conflicts,
        unique_features,
        alt_unique_features,
        overlap_graph,
        feature_list,
        alt_feature_list,
        file_name1,
        file_name2,
        outdir,
):
    """
    Write the summary and report files for each type of comparison.

    :param matching: List of 2-tuples of SeqFeatures of exactly matching annotations.
    :param conflicts: List of 2-tuples of SeqFeatures of conflicting annotations.
    :param unique_features: List of SeqFeatures containing unique features from the first annotation file.
    :param alt_unique_features: List of SeqFeatures containing unique features from the alternate annotation file.
    :param overlap_graph: networkx undirected graph with edges between overlapping features
    :param feature_list: List of features from the input.gbk file ordered by position.
    :param alt_feature_list: List of features from the alternative input.gbk file ordered by position.
    :param file_name1: String label for first annotation file.
    :param file_name2: String label for second annotation file.
    :param outdir: String name for the out directory
    """

    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except:
            sys.exit(f"Could not create directory {outdir}")

    path = f"{outdir}{os.sep}"
    summary_file = f"{path}summary.txt"
    pseudo_summary_file = f"{path}summary.pseudo.txt"
    matching_file = f"{path}colocated.tsv"
    conflicts_file = f"{path}conflicting.tsv"
    uniques_file = f"{path}{file_name1}.unique.tsv"
    alt_uniques_file = f"{path}{file_name2}.unique.tsv"

    matching_header = [
        "locus_tag_1", "gene_name_1", "pseudo_1", "pseudo_type1",
        "locus_tag_2", "gene_name_2", "pseudo_2", "pseudo_type2",
        "start", "end", "strand",
        "tags",
    ]
    conflicts_header = [
        "locus_tag_1", "gene_name_1", "start_1", "end_1", "strand_1", "pseudo_1", "pseudo_type1",
        "locus_tag_2", "gene_name_2", "start_2", "end_2", "strand_2", "pseudo_2", "pseudo_type2",
        "tags",
    ]
    uniques_header = [
        "locus_tag",
        "gene_name",
        "start",
        "end",
        "strand",
        "pseudo",
        "pseudo_type",
        "overlaps",
    ]
    files = [
        matching_file,
        conflicts_file,
        uniques_file, alt_uniques_file,
    ]
    headers = [
        matching_header,
        conflicts_header,
        uniques_header, uniques_header,
    ]
    reports = [
        matching,
        conflicts,
        unique_features, alt_unique_features,
    ]
    formatters = [
        format_match,
        format_conflict,
        format_unique, format_unique,
    ]

    for i in range(len(files)):
        with open(files[i], 'w') as f:
            print('\t'.join(headers[i]), file=f)
            for record in reports[i]:
                line = formatters[i](record, G=overlap_graph)
                print('\t'.join(line), file=f)

    pseudo_feature_list = [_ for _ in feature_list if designator.is_pseudo(_.qualifiers)]
    alt_pseudo_feature_list = [_ for _ in alt_feature_list if designator.is_pseudo(_.qualifiers)]
    pseudo_unique_features = [f for f in unique_features if designator.is_pseudo(f.qualifiers)]
    alt_pseudo_unique_features = [f for f in alt_unique_features if designator.is_pseudo(f.qualifiers)]
    pseudo_matching = [(f1, f2) for f1, f2 in matching if designator.is_pseudo(f1.qualifiers)]
    alt_pseudo_matching = [(f1, f2) for f1, f2 in matching if designator.is_pseudo(f2.qualifiers)]
    pseudo_conflicts = [(f1, f2) for f1, f2 in conflicts if designator.is_pseudo(f1.qualifiers)]
    alt_pseudo_conflicts = [(f1, f2) for f1, f2 in conflicts if designator.is_pseudo(f2.qualifiers)]

    breakdown = {
        'all': {
            'summary_file': summary_file,
            'file_name1': file_name1,
            'file_name2': file_name2,
            'total_features1': str(len(feature_list)),
            'total_features2': str(len(alt_feature_list)),
            'total_named_features1': str(len(get_named(feature_list))),
            'total_named_features2': str(len(get_named(alt_feature_list))),
            'unique1': str(len(unique_features)),
            'unique2': str(len(alt_unique_features)),
            'unique_named1': str(len(get_named(unique_features))),
            'unique_named2': str(len(get_named(alt_unique_features))),
            'matching1': str(len(matching)),
            'matching2': str(len(matching)),
            'matching_named1': str(len(get_named(matching, key=lambda _:_[0]))),
            'matching_named2': str(len(get_named(matching, key=lambda _:_[1]))),
            'matching_ident_named1': str(len(identically_named(matching))),
            'matching_ident_named2': str(len(identically_named(matching))),
            'matching_diff_named1': str(len(differentially_named(matching)[0])),
            'matching_diff_named2': str(len(differentially_named(matching)[1])),
            'conflicting_genes1': str(len(set([_[0].label for _ in conflicts]))),
            'conflicting_genes2': str(len(set([_[1].label for _ in conflicts]))),
            'conflicting_genes_named1': str(len(get_named(conflicts, key=lambda _:_[0], unique=True))),
            'conflicting_genes_named2': str(len(get_named(conflicts, key=lambda _:_[1], unique=True))),
            'conflicts1': str(len(conflicts)),
            'conflicts2': str(len(conflicts)),
            'conflicts_ident_named1': str(len(identically_named(conflicts))),
            'conflicts_ident_named2': str(len(identically_named(conflicts))),
            'conflicts_diff_named1': str(len(differentially_named(conflicts)[0])),
            'conflicts_diff_named2': str(len(differentially_named(conflicts)[1])),
        },
        'pseudo':{
            'summary_file': pseudo_summary_file,
            'file_name1': file_name1,
            'file_name2': file_name2,
            'total_features1': str(len(pseudo_feature_list)),
            'total_features2': str(len(alt_pseudo_feature_list)),
            'total_named_features1': str(len(get_named(pseudo_feature_list))),
            'total_named_features2': str(len(get_named(alt_pseudo_feature_list))),
            'unique1': str(len(pseudo_unique_features)),
            'unique2': str(len(alt_pseudo_unique_features)),
            'unique_named1': str(len(get_named(pseudo_unique_features))),
            'unique_named2': str(len(get_named(alt_pseudo_unique_features))),
            'matching1': str(len(pseudo_matching)),
            'matching2': str(len(alt_pseudo_matching)),
            'matching_named1': str(len(get_named(pseudo_matching, key=lambda _:_[0]))),
            'matching_named2': str(len(get_named(alt_pseudo_matching, key=lambda _:_[1]))),
            'matching_ident_named1': str(len(identically_named(pseudo_matching))),
            'matching_ident_named2': str(len(identically_named(alt_pseudo_matching))),
            'matching_diff_named1': str(len(differentially_named(pseudo_matching)[0])),
            'matching_diff_named2': str(len(differentially_named(alt_pseudo_matching)[1])),
            'conflicting_genes1': str(len(set([_[0].label for _ in pseudo_conflicts]))),
            'conflicting_genes2': str(len(set([_[1].label for _ in alt_pseudo_conflicts]))),
            'conflicting_genes_named1': str(len(get_named(pseudo_conflicts, key=lambda _:_[0], unique=True))),
            'conflicting_genes_named2': str(len(get_named(alt_pseudo_conflicts, key=lambda _:_[1], unique=True))),
            'conflicts1': str(len(pseudo_conflicts)),
            'conflicts2': str(len(alt_pseudo_conflicts)),
            'conflicts_ident_named1': str(len(identically_named(pseudo_conflicts))),
            'conflicts_ident_named2': str(len(identically_named(alt_pseudo_conflicts))),
            'conflicts_diff_named1': str(len(differentially_named(pseudo_conflicts)[0])),
            'conflicts_diff_named2': str(len(differentially_named(alt_pseudo_conflicts)[1])),
        }
    }

    for summary, results in breakdown.items():
        with open(results['summary_file'], 'w') as f:
            print('\t'.join([
                "",
                results['file_name1'],
                results['file_name2'],
            ]), file=f)
            print('\t'.join([
                f"total",
                results['total_features1'],
                results['total_features2'],
            ]), file=f)
            print('\t'.join([
                f"total_named",
                results['total_named_features1'],
                results['total_named_features2'],
            ]), file=f)
            #
            print(file=f)
            #
            print('\t'.join([
                f"unique",
                results['unique1'],
                results['unique2'],
            ]), file=f)
            print('\t'.join([
                f"unique_named",
                results['unique_named1'],
                results['unique_named2'],
            ]), file=f)
            #
            print(file=f)
            #
            print('\t'.join([
                f"colocated",
                results['matching1'],
                results['matching2'],
            ]), file=f)
            print('\t'.join([
                f"colocated_named",
                results['matching_named1'],
                results['matching_named2'],
            ]), file=f)
            print('\t'.join([
                f"colocated_identically_named",
                results['matching_ident_named1'],
                results['matching_ident_named2'],
            ]), file=f)
            print('\t'.join([
                f"colocated_exclusively_named",
                results['matching_diff_named1'],
                results['matching_diff_named2'],
            ]), file=f)
            #
            print(file=f)
            #
            print('\t'.join([
                f"conflicting_genes",
                results['conflicting_genes1'],
                results['conflicting_genes2'],
            ]), file=f)
            print('\t'.join([
                f"conflicting_genes_named",
                results['conflicting_genes_named1'],
                results['conflicting_genes_named2'],
            ]), file=f)
            print('\t'.join([
                f"conflicts",
                results['conflicts1'],
                results['conflicts2'],
            ]), file=f)
            print('\t'.join([
                f"conflicts_identically_named",
                results['conflicts_ident_named1'],
                results['conflicts_ident_named2'],
            ]), file=f)
            print('\t'.join([
                f"conflicts_exclusively_named",
                results['conflicts_diff_named1'],
                results['conflicts_diff_named2'],
            ]), file=f)
