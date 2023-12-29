from collections import defaultdict
import os
import sys

from intervaltree import Interval, IntervalTree
import networkx as nx

from .bio import SeqIO
from . import extractor
from . import fileManager
from . import annomerge
from . import designator


def main(args):

    gbk_file1, gbk_file2 = args.annotations
    outdir = args.outdir

    basename1, ext1 = os.path.splitext(os.path.basename(gbk_file1))
    basename2, ext2 = os.path.splitext(os.path.basename(gbk_file2))

    if (
            ext1[1:] not in fileManager.exts['genbank'] or
            ext2[1:] not in fileManager.exts['genbank']
    ):
        sys.exit(f"Input file must be in Genbank format")
    if basename1 == basename2:
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

    feature_list, pseudo_list = generate_record(gbk_file1)
    alt_feature_list, alt_pseudo_list = generate_record(gbk_file2)

    matching, conflicts, uniques, alt_uniques, ovl_graph = compare(
        feature_list,
        alt_feature_list,
    )

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

    pseudo_matching = [(f1, f2) for f1, f2 in matching
                       if designator.is_pseudo(f1.qualifiers) or designator.is_pseudo(f2.qualifiers)]
    pseudo_conflicts = [(f1, f2) for f1, f2 in conflicts
                       if designator.is_pseudo(f1.qualifiers) or designator.is_pseudo(f2.qualifiers)]
    pseudo_uniques = [f for f in uniques if designator.is_pseudo(f.qualifiers)]
    alt_pseudo_uniques = [f for f in alt_uniques if designator.is_pseudo(f.qualifiers)]

    write_reports(
        pseudo_matching,
        pseudo_conflicts,
        pseudo_uniques,
        alt_pseudo_uniques,
        ovl_graph,
        pseudo_list,
        alt_pseudo_list,
        basename1,
        basename2,
        outdir,
        "pseudo",
    )

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
    'ambiguous residues',
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

def generate_record(gbk):
    """
    Creates a defaultdictionary with all of the features from a .gbk annotation file.
    :param gbk: String of a path to a .gbk annotation file.
    :return anno_list: List of features from the input.gbk file ordered by position.
    :return pseudo_list: List of gene names with a 'pseudo' qualifier from the input.gbk file ordered by position.
    """
    anno_dict = defaultdict(list)
    pseudo_dict = defaultdict(list)
    anno_list = []
    pseudo_list = []
    all_records = SeqIO.parse(gbk, "genbank")

    for anno_record in all_records:
        np = which_np(anno_record)
        for f in anno_record.features:
            if f.type == 'CDS' or f.type == 'pseudo':
                if designator.is_pseudo(f.qualifiers):
                    f.evidence = np(f)
                    pseudo_list.append(f)
                anno_list.append(f)
    return anno_list, pseudo_list

def get_pseudo_type(feature):
    return f"{feature.evidence if hasattr(feature, 'evidence') else '.'}"

def compare(feature_list, alt_feature_list):
    """
    Organize each feature into various categories based on overlaps between the annotation files.

    :param feature_list: List of features from the input.gbk file ordered by position.
    :param alt_feature_list: List of features from the alternative input.gbk file ordered by position.
    :return matching: List of lists containing information for exactly matching annotations.
    :return conflicts: List of lists containing information for conflicting annotations.
    :return unique_features: List of lists containing information for unique features from the first annotation file.
    """

    def get_feature_interval(feature):
        return Interval(int(feature.location.start), int(feature.location.end))

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

    # A second bipartite graph to track all overlaps (except exact matches)
    # we don't want edge direction here because we want to easily identify edges originating from
    # features from either genome.
    G_all_partial_overlaps = G_conflict.copy().to_undirected()

    co_located = []
    conflicting = []
    unique = []
    alt_unique = []

    for i, feature in enumerate(feature_list):
        feature.label = f"X{i}"
        for G in [G_conflict, G_all_partial_overlaps]:
            G.add_node(
                feature.label,
                annotation=feature,
                bipartite=0,
            )
        overlaps = alt_interval_tree.overlap(get_feature_interval(feature))
        if overlaps:
            overlaps = [_[2] for _ in list(overlaps)]
            overlaps = sorted(overlaps, key=lambda _: _.location.start)

        colo = [_ for _ in overlaps if (annomerge.overlap_inframe(feature.location, _.location) and
                                               feature.location == _.location)]
        conf = [_ for _ in overlaps if (annomerge.overlap_inframe(feature.location, _.location) and
                                               feature.location != _.location)]
        non_conf = [_ for _ in overlaps if not annomerge.overlap_inframe(feature.location, _.location)]

        # 'co_located' : exact match
        if any(colo):
            # we don't care about further overlaps with these pairs, as they're fully accounted for.
            G_conflict.remove_nodes_from([_.label for _ in colo] + [feature.label])
            # the only time len(colo)>1 is if there are redundant annotation entries.
            _ = colo[0]
            co_located.append( (feature, colo[0]) )
            continue

        # 'conflicting' : overlaps in frame
        if any(conf):
            for alt_feature in conf:
                for G in [G_conflict, G_all_partial_overlaps]:
                    # alt_feature could have been dropped from the graph in a previous match.
                    # if so, don't add it back into the graph.
                    if alt_feature.label in G.nodes:
                        G.add_edge(feature.label, alt_feature.label, inframe=True)

        # 'non-conflicting' : no overlaps at all or overlaps out of frame
        if any(non_conf):
            G_all_partial_overlaps.add_edges_from(
                [(feature.label, _.label) for _ in non_conf],
                inframe=False,
            )

    #
    # See what's left now that the dust has settled
    #
    for comp in nx.weakly_connected_components(G_conflict):
        comp = list(comp)
        # "connected components" with only one node represent unique features.
        if len(comp) == 1:
            annotation = G_conflict.nodes[comp[0]]['annotation']
            if G_conflict.nodes[comp[0]]['bipartite'] == 0:
                unique.append(annotation)
            else:
                alt_unique.append(annotation)
        else:
            conflicting += [
                (G_conflict.nodes[x]['annotation'], G_conflict.nodes[y]['annotation'])
                for x,y in G_conflict.edges(comp)
            ]

    return co_located, conflicting, unique, alt_unique, G_all_partial_overlaps

#
# Helper functions for write_reports()
#

def unpack_feature(feature):
    return {
        'ltag': extractor.get_ltag(feature),
        'gene': extractor.get_gene(feature),
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
    if f1_pseudo and not f2_pseudo:
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

def get_named(feature_list, key=lambda _:_):
    """
    Get all features that have a name defined.

    :param feature_list: list of SeqFeatures or list of containers including SeqFeatures
    :param key: function that returns the SeqFeature from an item in feature_list
    """
    return [key(_) for _ in feature_list if name(key(_))]

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
        suffix="",
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
    :param suffix: Optional argument to change the suffix of report files.
    """

    features_total = len(feature_list)
    alt_features_total = len(alt_feature_list)

    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except:
            sys.exit(f"Could not create directory {outdir}")

    path = f"{outdir}{os.sep}"
    summary_file = f"{path}summary{'.' if suffix else ''}{suffix}.txt"
    matching_file = f"{path}colocated{'.' if suffix else ''}{suffix}.tsv"
    conflicts_file = f"{path}conflicting{'.' if suffix else ''}{suffix}.tsv"
    uniques_file = f"{path}{file_name1}.unique{'.' if suffix else ''}{suffix}.tsv"
    alt_uniques_file = f"{path}{file_name2}.unique{'.' if suffix else ''}{suffix}.tsv"

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

    colo_same_names = identically_named(matching)
    colo_exc_named, alt_colo_exc_named = differentially_named(matching)

    confl_same_names = identically_named(conflicts)
    confl_exc_named, alt_confl_exc_named = differentially_named(conflicts)

    for i in range(len(files)):
        with open(files[i], 'w') as f:
            print('\t'.join(headers[i]), file=f)
            for record in reports[i]:
                line = formatters[i](record, G=overlap_graph)
                print('\t'.join(line), file=f)

    uniq_c = len(set([_[0].label for _ in conflicts]))
    alt_uniq_c = len(set([_[1].label for _ in conflicts]))

    with open(summary_file, 'w') as f:
        print('\t'.join(["", file_name1, file_name2]), file=f)
        print('\t'.join([f"total", str(features_total), str(alt_features_total)]), file=f)
        print('\t'.join([
            f"total_named",
            str(len(get_named(feature_list))),
            str(len(get_named(alt_feature_list))),
        ]), file=f)
        #
        print(file=f)
        #
        print('\t'.join([f"unique", str(len(unique_features)), str(len(alt_unique_features))]), file=f)
        print('\t'.join([
            f"unique_named",
            str(len(get_named(unique_features))),
            str(len(get_named(alt_unique_features))),
        ]), file=f)
        #
        print(file=f)
        #
        print('\t'.join([f"co-located", str(len(matching)), str(len(matching))]), file=f)
        print('\t'.join([
            f"co-located_named",
            str(len(get_named(matching, key=lambda _:_[0]))),
            str(len(get_named(matching, key=lambda _:_[1]))),
        ]), file=f)
        print('\t'.join([
            f"co-located_identically_named",
            str(len(colo_same_names)),
            str(len(colo_same_names)),
        ]), file=f)
        print('\t'.join([
            f"co-located_exclusively_named",
            str(len(colo_exc_named)),
            str(len(alt_colo_exc_named)),
        ]), file=f)
        #
        print(file=f)
        #
        print('\t'.join([f"conflicting", str(uniq_c), str(alt_uniq_c)]), file=f)
        print('\t'.join([
            f"conflicting_named",
            str(len(get_named(conflicts, key=lambda _:_[0]))),
            str(len(get_named(conflicts, key=lambda _:_[1]))),
        ]), file=f)
        print('\t'.join([
            f"conflicting_identically_named",
            str(len(confl_same_names)),
            str(len(confl_same_names)),
        ]), file=f)
        print('\t'.join([
            f"conflicting_exclusively_named",
            str(len(confl_exc_named)),
            str(len(alt_confl_exc_named)),
        ]), file=f)
