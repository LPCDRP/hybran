import os
import csv
import sys
from collections import defaultdict
from . import extractor
from . import fileManager
from . import annomerge
from . import designator
from .bio import SeqIO

def main(args):

    gbk_file1 = args.annotations[0]
    gbk_file2 = args.annotations[1]
    outdir = args.outdir

    feature_list, pseudo_list = generate_record(gbk_file1)
    alt_feature_list, alt_pseudo_list = generate_record(gbk_file2)

    matching, conflicts, uniques = compare(feature_list, alt_feature_list)
    alt_matching, alt_conflicts, alt_uniques = compare(alt_feature_list, feature_list)

    pseudo_matching, pseudo_conflicts, pseudo_uniques = compare(pseudo_list, alt_feature_list)
    alt_pseudo_matching, alt_pseudo_conflicts, alt_pseudo_uniques = compare(alt_pseudo_list, feature_list)

    write_reports(matching, alt_matching, conflicts, alt_conflicts, uniques, alt_uniques,
                  feature_list, alt_feature_list, gbk_file1, gbk_file2, outdir)

    write_reports(pseudo_matching, alt_pseudo_matching, pseudo_conflicts, alt_pseudo_conflicts, pseudo_uniques, alt_pseudo_uniques,
                  pseudo_list, alt_pseudo_list, gbk_file1, gbk_file2, outdir, "pseudo")

def furthest_location(loc1, loc2):
    """
    Finds the SeqFeature.Location that is the farthest downstream.
    :param loc1: SeqFeature location.
    :param loc2: SeqFeature location.
    :return: Downstream SeqFeature location or None if the locations are exactly the same.
    """
    same_start = False
    same_stop = False
    if loc1.start == loc2.start:
        same_start = True
    if loc1.end == loc2.end:
        same_stop = True

    if same_start and same_stop:
        return None
    elif not same_start and same_stop:
        if loc1.start > loc2.start:
            return loc1
        return loc2
    else:
        if loc1.end > loc2.end:
            return loc1
        return loc2

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
    print(f"Generating record for {gbk.split('/')[-1]}")

    for anno_record in all_records:
        np = which_np(anno_record)
        for f in anno_record.features:
            if f.type == 'CDS' or f.type == 'pseudo':
                gene_name = (f.qualifiers.get('locus_tag')[0] if f.qualifiers.get('gene') == None else f.qualifiers.get('gene')[0])
                #Known bug from RATT annotations creates a nonexistent copy of mamB in hybran
                if gene_name == 'mamB' and "RATT" in f.qualifiers['inference'][0]:
                    continue

                if designator.is_pseudo(f.qualifiers):
                    f.record = np(f)
                    pseudo_list.append(f)
                anno_list.append(f)

    return anno_list, pseudo_list

def is_overlap(loc1, loc2):
    """
    Determine if two locations overlap.
    :param loc1: SeqFeature location.
    :param loc2: SeqFeature location.
    :return: True if two SeqFeature locations overlap in any way (in frame or out of frame). Otherwise return False.
    """
    overlap = (min(loc1.end, loc2.end) - max(loc1.start, loc2.start))
    if loc1.strand == loc2.strand and overlap > 0:
        return True
    return False

def compare(feature_list, alt_feature_list):
    """
    Organize each feature into various categories based on overlaps between the annotation files.

    :param feature_list: List of features from the input.gbk file ordered by position.
    :param alt_feature_list: List of features from the alternative input.gbk file ordered by position.
    :return matching: List of lists containing information for exactly matching annotations.
    :return conflicts: List of lists containing information for conflicting annotations.
    :return unique_features: List of lists containing information for unique features from the first annotation file.
    """
    def get_pseudo_type(feature):
        return f"{feature.record if hasattr(feature, 'record') else '.'}"
    
    def get_surrounding_features(feature, alt_feature_list):
        """
        :param feature: A SeqFeature from the first annotation file
        :param alt_feature_list: List of SeqFeatures from the alternate annotation file
        :return overlaps: A list of upstream and downstream features from the alternate annotation file
        that overlap with the input feature
        """
        overlaps = []
        next_index = 0

        for i in range(len(alt_feature_list)):
            if is_overlap(feature.location, alt_feature_list[i].location):
                overlaps.append(alt_feature_list[i])
        return overlaps

    #Organize each feature into 3 categories: 'matching', 'unique', 'conflicts'
    #Finds uniques from the first annotation file.
    matching = []
    conflicts = []
    unique_features = []
    for feature in feature_list:
        frame = []
        locus_tag = extractor.get_ltag(feature)
        gene_name = extractor.get_gene(feature)
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        pseudo = int(designator.is_pseudo(feature.qualifiers))
        overlaps = get_surrounding_features(feature, alt_feature_list)

        if len(overlaps) > 0:
            for i in range(len(overlaps)):
                frame.append(annomerge.overlap_inframe(feature.location, overlaps[i].location))

        break_check = False
        for i in range(len(frame)):
            # 'matching' : exact match
            if frame[i] == True and (str(feature.location) == str(overlaps[i].location)):
                line = [
                    locus_tag, gene_name, pseudo, get_pseudo_type(feature),
                    extractor.get_ltag(overlaps[i]),
                    extractor.get_gene(overlaps[i]),
                    int(designator.is_pseudo(overlaps[i].qualifiers)),
                    get_pseudo_type(overlaps[i]),
                    start, end, strand,
                ]
                matching.append(line)
                break_check = True
                break
        if break_check:
            continue

        for i in range(len(frame)):
            # 'conflicting' : overlaps in frame
            if frame[i] == True and (str(feature.location) != str(overlaps[i].location)):
                line = [
                    locus_tag, gene_name, start, end, strand, pseudo, get_pseudo_type(feature),
                    extractor.get_ltag(overlaps[i]),
                    extractor.get_gene(overlaps[i]),
                    int(overlaps[i].location.start),
                    int(overlaps[i].location.end),
                    overlaps[i].location.strand,
                    int(designator.is_pseudo(overlaps[i].qualifiers)),
                    get_pseudo_type(overlaps[i])
                ]
                conflicts.append(line)
                break_check = True
        if break_check:
            continue

        # 'unique' : no overlaps or overlaps out of frame
        if not any(overlaps) or False in frame:
            overlap_note = [
                f"{extractor.get_gene(_)}|{str(_.location)}|pseudo={int(designator.is_pseudo(_.qualifiers))}|"
                f"pseudo_type={get_pseudo_type(_)}" for _ in overlaps]
            line = [locus_tag, gene_name, start, end, strand, pseudo, get_pseudo_type(feature), overlap_note]
            unique_features.append(line)

    print(f"Found all matches, conflicts, and unique features.")
    return matching, conflicts, unique_features

def write_reports(matching, alt_matching, conflicts, alt_conflicts, unique_features, alt_unique_features,
                  feature_list, alt_feature_list, gbk_file1, gbk_file2, outdir, suffix=""):
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
    :param gbk_file1: String of a path to a .gbk annotation file.
    :param gbk_file2: String of a path to an alternate .gbk annotation file.
    :param outdir: String name for the out directory
    :param suffix: Optional argument to change the suffix of report files.
    """

    file_name1 = os.path.splitext(os.path.basename(gbk_file1))[0]
    file_name2 = os.path.splitext(os.path.basename(gbk_file2))[0]
    features_total = len(feature_list)
    alt_features_total = len(alt_feature_list)

    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except:
            sys.exit(f"Could not create directory {outdir}")

    summary_file = f"{outdir}/summary{'_' if suffix else ''}{suffix}.txt"

    matching_file = f"{outdir}/matching{'_' if suffix else ''}{suffix}.csv"
    matching_header = [
                "locus_tag_1", "gene_name_1", "pseudo_1", "pseudo_type1",
                "locus_tag_2", "gene_name_2", "pseudo_2", "pseudo_type2",
                "start", "end", "strand",
    ]
    conflicts_file = f"{outdir}/conflicts{'_' if suffix else ''}{suffix}.csv"
    conflicts_header = [
                "locus_tag_1", "gene_name_1", "start_1", "end_1", "strand_1", "pseudo_1", "pseudo_type1",
                "locus_tag_2", "gene_name_2", "start_2", "end_2", "strand_2", "pseudo_2", "pseudo_type2",
    ]
    uniques_file = f"{outdir}/{file_name1}_uniques{'_' if suffix else ''}{suffix}.csv"
    alt_uniques_file = f"{outdir}/{file_name2}_uniques{'_' if suffix else ''}{suffix}.csv"
    uniques_header = [
        "locus_tag", "gene_name", "start",
        "end", "strand", "pseudo", "pseudo_type", "overlaps"
    ]


    files = [matching_file, conflicts_file, uniques_file, alt_uniques_file]
    headers = [matching_header, conflicts_header, uniques_header, uniques_header]
    reports = [matching, conflicts, unique_features, alt_unique_features]

    for i in range(len(files)):
        with open(files[i], 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            header = headers[i]
            writer.writerow(header)
            for line in reports[i]:
                writer.writerow(line)

    uniq_f = str(len(conflicts) - len(set([str(_[0:6]) for _ in conflicts])))
    alt_uniq_f = str(len(alt_conflicts) - len(set([str(_[0:6]) for _ in alt_conflicts])))

    with open(summary_file, 'w') as f:
        print('\t'.join([f"", f"{file_name1}", f"{file_name2}"]), file=f)
        print('\t'.join([f"Total Features", str(features_total), str(alt_features_total)]), file=f)
        print('\t'.join([f"Number of Unique Features", str(len(unique_features)), str(len(alt_unique_features))]), file=f)
        print('\t'.join([f"Number of Exact Matches", str(len(matching)), str(len(alt_matching))]), file=f)
        print('\t'.join([f"Number of Conflicts", str(len(conflicts)), str(len(alt_conflicts))]), file=f)
        print('\t'.join([f"Unique Conflicts", uniq_f, alt_uniq_f]), file=f)
