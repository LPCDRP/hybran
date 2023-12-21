from collections import defaultdict
import os
import sys

from intervaltree import Interval, IntervalTree

from .bio import SeqIO
from . import extractor
from . import fileManager
from . import annomerge
from . import designator


def main(args):

    gbk_file1, gbk_file2 = args.annotations
    outdir = args.outdir

    basename1 = os.path.splitext(os.path.basename(gbk_file1))[0]
    basename2 = os.path.splitext(os.path.basename(gbk_file2))[0]

    if (
            os.path.splitext(os.path.basename(gbk_file1))[1] != ".gbk" or
            os.path.splitext(os.path.basename(gbk_file2))[1] != ".gbk"
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

    matching, conflicts, uniques = compare(feature_list, alt_feature_list)
    alt_matching, alt_conflicts, alt_uniques = compare(alt_feature_list, feature_list)

    pseudo_matching, pseudo_conflicts, pseudo_uniques = compare(pseudo_list, alt_feature_list)
    alt_pseudo_matching, alt_pseudo_conflicts, alt_pseudo_uniques = compare(alt_pseudo_list, feature_list)

    write_reports(
        matching,
        alt_matching,
        conflicts,
        alt_conflicts,
        uniques,
        alt_uniques,
        feature_list,
        alt_feature_list,
        basename1,
        basename2,
        outdir,
    )

    write_reports(
        pseudo_matching,
        alt_pseudo_matching,
        pseudo_conflicts,
        alt_pseudo_conflicts,
        pseudo_uniques,
        alt_pseudo_uniques,
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
    evidence_code = pseudo_notes[0].split(":", 2)[2]
    return evidence_code

def pgap_np(feature):
    """
    Note parser to identify types of pseudogenes specifically from PGAP annotation files
    :param feature: A SeqFeature
    :return evidence_code: String from the note field of a CDS with a '/pseudo' qualifier
    """
    evid = {}
    evidence_code = []
    evid['frameshift'] = any([_ for _ in feature.qualifiers['note'] if "frameshifted" in _])
    evid['internal_stop'] = any([_ for _ in feature.qualifiers['note'] if "internal_stop" in _])
    evid['incomplete'] = any([_ for _ in feature.qualifiers['note'] if "incomplete" in _])
    for k, v in evid.items():
        if v:
            evidence_code.append(str(k))
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

    def get_pseudo_type(feature):
        return f"{feature.evidence if hasattr(feature, 'evidence') else '.'}"

    alt_interval_tree = IntervalTree()
    for alt_feature in alt_feature_list:
        alt_interval_tree.addi(int(alt_feature.location.start), int(alt_feature.location.end), alt_feature)

    #Organize each feature into 3 categories: 'co_located', 'conflicting', 'non_conflicting'
    #Finds uniques from the first annotation file.
    co_located = []
    conflicting = []
    non_conflicting = []

    for feature in feature_list:
        if feature.type != 'CDS':
            continue
        locus_tag = extractor.get_ltag(feature)
        gene_name = extractor.get_gene(feature)
        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand
        pseudo = int(designator.is_pseudo(feature.qualifiers))
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
            _ = colo[0]
            line = [
                locus_tag, gene_name, pseudo, get_pseudo_type(feature),
                extractor.get_ltag(_), extractor.get_gene(_),
                int(designator.is_pseudo(_.qualifiers)),
                get_pseudo_type(_), start, end, strand,
            ]
            co_located.append(line)
            continue

        # 'conflicting' : overlaps in frame
        if any(conf):
            for _ in conf:
                line = [
                    locus_tag, gene_name, start, end, strand, pseudo, get_pseudo_type(feature),
                    extractor.get_ltag(_), extractor.get_gene(_), (_.location.start),
                    (_.location.end), _.location.strand, int(designator.is_pseudo(_.qualifiers)),
                    get_pseudo_type(_)
                ]
                conflicting.append(line)
            continue

        # 'non-conflicting' : no overlaps at all or overlaps out of frame
        overlap_note = "."
        if any(non_conf):
            overlap_note = ";".join([
                f"{extractor.get_gene(_)}|{str(_.location)}|pseudo={int(designator.is_pseudo(_.qualifiers))}|"
                f"pseudo_type={get_pseudo_type(_)}" for _ in non_conf
            ])
        line = [locus_tag, gene_name, start, end, strand, pseudo, get_pseudo_type(feature), overlap_note]
        non_conflicting.append(line)
    return co_located, conflicting, non_conflicting

def write_reports(
        matching,
        alt_matching,
        conflicts,
        alt_conflicts,
        unique_features,
        alt_unique_features,
        feature_list,
        alt_feature_list,
        file_name1,
        file_name2,
        outdir,
        suffix="",
):
    """
    Write the summary and report files for each type of comparison.

    :param matching: List of lists containing information for exactly matching annotations.
    :param alt_matching: List of lists containing information for exactly matching annotations.
    :param conflicts: List of lists containing information for conflicting annotations from the first annotation file.
    :param alt_conflicts: List of lists containing information for conflicting annotations from the alternate annotation file.
    :param unique_features: List of lists containing information for unique features from the first annotation file.
    :param alt_unique_features: List of lists containing information for unique features from the alternate annotation file.
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

    path = f"{outdir}/"
    summary_file = f"{path}summary{'.' if suffix else ''}{suffix}.txt"
    matching_file = f"{path}{file_name1}.colo{'.' if suffix else ''}{suffix}.tsv"
    alt_matching_file = f"{path}{file_name2}.colo{'.' if suffix else ''}{suffix}.tsv"
    conflicts_file = f"{path}{file_name1}.confl{'.' if suffix else ''}{suffix}.tsv"
    alt_conflicts_file = f"{path}{file_name2}.confl{'.' if suffix else ''}{suffix}.tsv"
    uniques_file = f"{path}{file_name1}.nonconfl{'.' if suffix else ''}{suffix}.tsv"
    alt_uniques_file = f"{path}{file_name2}.nonconfl{'.' if suffix else ''}{suffix}.tsv"

    matching_header = [
        "locus_tag_1", "gene_name_1", "pseudo_1", "pseudo_type1",
        "locus_tag_2", "gene_name_2", "pseudo_2", "pseudo_type2",
        "start", "end", "strand",
    ]
    conflicts_header = [
        "locus_tag_1", "gene_name_1", "start_1", "end_1", "strand_1", "pseudo_1", "pseudo_type1",
        "locus_tag_2", "gene_name_2", "start_2", "end_2", "strand_2", "pseudo_2", "pseudo_type2",
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
        matching_file, alt_matching_file,
        conflicts_file, alt_conflicts_file,
        uniques_file, alt_uniques_file,
    ]
    headers = [
        matching_header, matching_header,
        conflicts_header, conflicts_header,
        uniques_header, uniques_header,
    ]
    reports = [
        matching, alt_matching,
        conflicts, alt_conflicts,
        unique_features, alt_unique_features,
    ]

    for i in range(len(files)):
        with open(files[i], 'w') as f:
            print('\t'.join(headers[i]), file=f)
            for line in reports[i]:
                line = [str(_) for _ in line]
                print('\t'.join(line), file=f)

    uniq_c = len(set([tuple(_[0:6]) for _ in conflicts]))
    alt_uniq_c = len(set([tuple(_[0:6]) for _ in alt_conflicts]))

    with open(summary_file, 'w') as f:
        print('\t'.join(["", file_name1, file_name2]), file=f)
        print('\t'.join([f"Total", str(features_total), str(alt_features_total)]), file=f)
        print('\t'.join([f"Non-conflicting", str(len(unique_features)), str(len(alt_unique_features))]), file=f)
        print('\t'.join([f"Co-located", str(len(matching)), str(len(alt_matching))]), file=f)
        print('\t'.join([f"Conflicting", str(uniq_c), str(alt_uniq_c)]), file=f)
