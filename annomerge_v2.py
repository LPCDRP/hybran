#!/usr/bin/env python2.7

__author__ = "Deepika Gunasekaran"
__version__ = "1.0.1"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@sdsu.edu"
__status__ = "Development"

# Title: Merge annotation from Prokka for the positions which are not annotated by RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files in case of
# multiple contigs/chromosome annotations and a Genbank file (.gbf) file from Prokka annotation run with a reference and
# a Genbank file (.gbf) file from Prokka annotation run without a reference. The output is an EMBL file with annotation
# predominantly from RATT and the intergenic regions annotated by RATT are filled with Prokka. This script also
# generates a log file to indicate characteristics of the transferred features from Prokka.

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import translate
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SeqRecord import SeqRecord
import collections
from numpy import median
import os
import pickle
import datetime


def get_ordered_features(feature_list):
    """
    :param feature_list: list of features of type SeqFeature
    :return: sorted list of features, sorted by genomic location
    """

    features_dict = {}
    ordered_features = []
    for feature in feature_list:
        feature_start = int(feature.location.start)
        if feature_start not in features_dict.keys():
            features_dict[feature_start] = [feature]
        else:
            features_dict[feature_start].append(feature)
    ordered_features_dict = collections.OrderedDict(sorted(features_dict.items()))
    for position in ordered_features_dict.keys():
        features_in_position = ordered_features_dict[position]
        for feature in features_in_position:
            ordered_features.append(feature)
    return ordered_features


def remove_duplicate_cds(feature_list):
    """
    :param feature_list:  list of features of type SeqFeature
    :return: list of features after removing duplicate CDSs i.e. annotation of same gene in same position
    """
    unique_feature_list = []
    added_cds = []
    for feature in feature_list:
        if feature.type != 'CDS':
           unique_feature_list.append(feature)
        else:
            feature_key = feature.qualifiers['locus_tag'][0] + ':' + str(int(feature.location.start)) + ':' + \
                          str(int(feature.location.end)) + ':' + str(int(feature.location.strand))
            if feature_key in added_cds:
                continue
            else:
                added_cds.append(feature_key)
                unique_feature_list.append(feature)
    return unique_feature_list


def generate_feature_dictionary(feature_list):
    # This function takes as input a list of features and returns a dictionary with the key as a tuple of feature
    # start and stop positions and the value as the feature.
    """
    :param feature_list: List of features (SeqFeature objects)
    :return: sorted dictionary ordered by the genomic position i.e. feature location
    """
    feature_dict = {}
    for feature in feature_list:
        feature_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
        feature_dict[feature_key] = feature
    sorted_feature_dict = collections.OrderedDict(sorted(feature_dict.items()))
    return sorted_feature_dict


def get_prokka_features_with_gene_name(all_prokka_features):
    """
    :param all_prokka_features: Input is a list of features from Prokka
    :return: returns 2 dictionarys; prokka annotations corresponding to Rv, and novel annotations where keys are locus
    tags and values are features
    """
    prokka_rv_features = {}
    novel_cds_prokka = {}
    for feature in all_prokka_features:
        if feature.type == 'CDS' and 'gene' in feature.qualifiers.keys():
            rv_gene = feature.qualifiers['gene'][0]
            if rv_gene not in prokka_rv_features.keys():
                prokka_rv_features[rv_gene] = [feature]
            else:
                prokka_rv_features[rv_gene].append(feature)
        elif feature.type == 'CDS' and 'gene' not in feature.qualifiers.keys():
            locus_tag = feature.qualifiers['locus_tag'][0]
            if locus_tag not in novel_cds_prokka.keys():
                novel_cds_prokka[locus_tag] = [feature]
            else:
                novel_cds_prokka[locus_tag].append(feature)
    output_file.write('Prokka annotation from Rv: ' + str(len(prokka_rv_features.keys())) + '\n')
    output_file.write('Prokka novel annotations: ' + str(len(novel_cds_prokka)) + '\n')
    return prokka_rv_features, novel_cds_prokka


def check_for_dnaA(feature_list):
    # This function takes as input a list of features and checks if the genome was circularized
    """
    :param feature_list: List of features in the RATT annotation record (EMBL file)
    :return: None
    """
    feature_dictionary = generate_feature_dictionary(feature_list)
    for cds in feature_dictionary.keys():
        if feature_dictionary[cds].type != 'CDS':
            continue
        if feature_dictionary[cds].qualifiers['locus_tag'][0] != 'Rv0001':
            print('DnaA is not the first gene in this genome. Possible circularization error')
            sys.exit()
        elif feature_dictionary[cds].qualifiers['locus_tag'][0] == 'Rv0001' and \
                int(feature_dictionary[cds].location.start) == 0:
            break
        else:
            print('DnaA is the first gene in this genome. But the position is off')
            sys.exit()
    return


def get_ratt_corrected_genes(ratt_report_fp):
    """
    :param ratt_report_fp: Correction Report from RATT
    :return: returns list of genes whose start/stop has been corrected by RATT
    """
    corrected_genes_list = []
    report_raw = open(ratt_report_fp).readlines()
    for line in report_raw:
        if len(line) <= 1:
            continue
        gene = line.strip().split()[0]
        if gene[:2] != 'Rv':
            continue
        if gene not in corrected_genes_list:
            corrected_genes_list.append(gene)
    return corrected_genes_list


def get_sorted_strand_genes(feature_list, strand):
    """
    :param feature_list: List of features (SeqFeature objects)
    :param strand: genomic strand for which features are to be sorted
    :return: dictionary of sorted CDSs in the starnd where the key is start of the CDS and value is a tuple of
    (end of CDS, gene name)
    """
    unsorted_features = {}
    for feature in feature_list:
        if int(feature.location.strand) != strand:
            continue
        elif feature.type != 'CDS':
            continue
        else:
            cds_start = int(feature.location.start)
            cds_end = int(feature.location.end)
            gene_name = feature.qualifiers['locus_tag'][0]
            unsorted_features[cds_start] = (cds_end, gene_name)
    sorted_features = collections.OrderedDict(sorted(unsorted_features.items()))
    return sorted_features


def rename_locus(gene, strand):
    """
    :param gene: gene name to be modified
    :param strand: strand information to check if gene is present in complementary strand
    :return: new gene name for merged gene
    """
    # This function checks names of existing locus tags in H37Rv and names the newly merged gene ensuring that the new
    # name does not conflict with existing locus tags
    char_start = 65 # ASCII conversion of 'A'
    gene_name = gene[:6]
    if strand == '-1' or strand == '-' or strand == -1:
        # Adding in 'c' to denote gene is located in the complementary strand
        new_gene_name = gene_name + chr(char_start) + 'c'
    else:
        new_gene_name = gene_name + chr(char_start)
    while new_gene_name in genes_in_rv:
        # If the assigned new gene name is in H37Rv, increment the alphabet suffix
        char_start += 1
        if strand == '-1' or strand == '-' or strand == -1:
            new_gene_name = gene_name + chr(char_start) + 'c'
        else:
            new_gene_name = gene_name + chr(char_start)
    return new_gene_name


def identify_merged_genes(ratt_features):
    # This function takes as input a list features annotated in RATT and identifies instances of merged CDS annotations.
    # The function returns a Boolean value (True if such instances are identified) and a dictionary with the strand as
    # the key and locations of the merged genes as values
    ratt_annotations = {}   # This dictionary holds 2 keys '-1' and '+1' denoting the strand and locations of all CDS
    # annotations in RATT for each strand
    ratt_unmerged_genes = {}
    ratt_merged_genes = {}
    merged_genes = False
    if len(ratt_features) == 0:
        return [], merged_genes
    for feature in ratt_features:
        # These are the features with 'joins'
        if 'Bio.SeqFeature.CompoundLocation' in str(type(feature.location)):
            continue
        # Getting all locations of CDS annotations and storing them in ratt_annotations
        if feature.location.strand not in ratt_annotations.keys() and feature.type == 'CDS':
            ratt_annotations[feature.location.strand] = [(int(feature.location.start), int(feature.location.end))]
            ratt_merged_genes[feature.location.strand] = []
            ratt_unmerged_genes[feature.location.strand] = []
        elif feature.type == 'CDS':
            ratt_annotations[feature.location.strand].append((int(feature.location.start), int(feature.location.end)))
    for strand in ratt_annotations.keys():
        for gene_location in ratt_annotations[strand]:
            if gene_location not in ratt_unmerged_genes[strand] and gene_location not in ratt_merged_genes[strand]:
                # First instance of the CDS location
                ratt_unmerged_genes[strand].append(gene_location)
            elif gene_location in ratt_unmerged_genes[strand] and gene_location not in ratt_merged_genes[strand]:
                # Identified duplicate CDS location annotation
                ratt_merged_genes[strand].append(gene_location)
                ratt_unmerged_genes[strand].remove(gene_location)
                merged_genes = True
            elif gene_location not in ratt_unmerged_genes[strand] and gene_location in ratt_merged_genes[strand]:
                # To account for instances where more than 1 gene is merged
                merged_genes = True
                continue
            else:
                print('UNACCOUNTED CASE')
                print(gene_location)
    return ratt_merged_genes, merged_genes


def get_annotation_for_merged_genes(merged_genes, prokka_features, ratt_features):
    # This function takes as input the dictionary of merged genes from the identify_merged_genes function, a list of
    # Prokka features and a list of RATT features and returns annotation for the merged genes from,
    # Prokka: If the genes are annotated as merged in both RATT and Prokka
    # RATT: If Prokka does not annotate these genes as merged
    check_positions_for_annotation = dict(merged_genes)
    merged_gene_locus = {}
    merged_features_addition = []
    features_from_ratt = {}
    genes_from_ratt = {}
    final_ratt_features = []
    final_prokka_features = []
    locus_to_remove_gene_tags = []
    for feature in ratt_features:
        # Identify CDS tags that should not be included in final_ratt_features
        if feature.type != 'CDS':
            continue
        ratt_feature_location = (int(feature.location.start), int(feature.location.end))
        if ratt_feature_location not in check_positions_for_annotation[feature.strand]:
            final_ratt_features.append(feature)
        elif ratt_feature_location in check_positions_for_annotation[feature.strand] and \
                        ratt_feature_location not in merged_gene_locus.keys():
            merged_gene_locus[ratt_feature_location] = [feature.qualifiers['locus_tag'][0]]
            locus_to_remove_gene_tags.append(feature.qualifiers['locus_tag'][0])
        elif ratt_feature_location in check_positions_for_annotation[feature.strand] and \
                        ratt_feature_location in merged_gene_locus.keys():
            merged_gene_locus[ratt_feature_location].append(feature.qualifiers['locus_tag'][0])
            locus_to_remove_gene_tags.append(feature.qualifiers['locus_tag'][0])
    for other_features in ratt_features:
        # Identifies 'gene' tags which corresponds to merged CDSs and does not include them in the final annotation
        if other_features.type == 'CDS':
            continue
        elif other_features.type == 'gene':
            if 'locus_tag' in other_features.qualifiers.keys():
                if other_features.qualifiers['locus_tag'][0] not in locus_to_remove_gene_tags:
                    final_ratt_features.append(other_features)
            else:
                print('This gene feature does not have a locus tag annotation in RATT')
                print(other_features)
        else:
            final_ratt_features.append(other_features)
    if len(prokka_features) == 0:
        for feature in ratt_features:
            if len(merged_genes[feature.location.strand]) == 0:
                continue
            merged_genes_in_strand = merged_genes[feature.location.strand]
            feature_location = (int(feature.location.start), int(feature.location.end))
            if feature_location in merged_genes_in_strand and feature.type == 'CDS':
                if feature_location not in features_from_ratt.keys():
                    features_from_ratt[feature_location] = [feature]
                    genes_from_ratt[feature_location] = [feature.qualifiers['locus_tag'][0]]
                elif feature.type == 'CDS':
                    features_from_ratt[feature_location].append(feature)
                    genes_from_ratt[feature_location].append(feature.qualifiers['locus_tag'][0])
    else:
        for feature in prokka_features:
            if len(merged_genes[feature.location.strand]) == 0:
                final_prokka_features.append(feature)
                continue
            merged_genes_in_strand = merged_genes[feature.location.strand]
            feature_location = (int(feature.location.start), int(feature.location.end))
            if feature_location in merged_genes_in_strand and feature.type == 'CDS':
                merged_features_string = ",".join(merged_gene_locus[feature_location])
                if 'note' not in feature.qualifiers.keys():
                    feature.qualifiers['note'] = ['The genes ' + merged_features_string +
                                                  ' in H37Rv(NC_000962.3) are merged in this isolate '
                                                  '(annotation from Prokka)']
                else:
                    feature.qualifiers['note'].append('The genes ' + merged_features_string +
                                                      ' in H37Rv(NC_000962.3) are merged in this isolate '
                                                      '(annotation from Prokka)')
                merged_genes[feature.location.strand].remove(feature_location)
                merged_features_addition.append(feature)
            else:
                final_prokka_features.append(feature)
        for feature in ratt_features:
            if len(merged_genes[feature.location.strand]) == 0:
                continue
            merged_genes_in_strand = merged_genes[feature.location.strand]
            feature_location = (int(feature.location.start), int(feature.location.end))
            if feature_location in merged_genes_in_strand and feature.type == 'CDS':
                if feature_location not in features_from_ratt.keys():
                    features_from_ratt[feature_location] = [feature]
                    genes_from_ratt[feature_location] = [feature.qualifiers['locus_tag'][0]]
                elif feature.type == 'CDS':
                    features_from_ratt[feature_location].append(feature)
                    genes_from_ratt[feature_location].append(feature.qualifiers['locus_tag'][0])
    for location in features_from_ratt.keys():
        new_feature = features_from_ratt[location][0]
        merged_features_string = ",".join(genes_from_ratt[location])
        new_gene_name = rename_locus(new_feature.qualifiers['locus_tag'][0], new_feature.location.strand)
        new_feature.qualifiers['locus_tag'] = [new_gene_name]
        if 'note' in new_feature.qualifiers.keys():
            new_feature.qualifiers['note'].append('The genes ' + merged_features_string +
                                                  ' in H37Rv(NC_000962.3) are merged in this isolate '
                                                  '(annotation from RATT)')
        else:
            new_feature.qualifiers['note'] = ['The genes ' + merged_features_string +
                                              ' in H37Rv(NC_000962.3) are merged in this isolate '
                                              '(annotation from RATT)']
        merged_features_addition.append(new_feature)
        merged_genes[new_feature.location.strand].remove(location)
    if len(merged_genes[-1]) > 0 or len(merged_genes[1]) > 0:
        corner_cases = True
    else:
        corner_cases = False
    return merged_features_addition, corner_cases, merged_genes, final_ratt_features, final_prokka_features


def check_corresponding_prokka_annotation(overlapping_gene_info, prokka_cds):
    locus_tag = overlapping_gene_info[0]
    ratt_gene_start = overlapping_gene_info[1]
    ratt_gene_end = overlapping_gene_info[2]
    ratt_gene_strand = overlapping_gene_info[3]
    return_feature_dict = {}
    prokka_cds_keys = '\t'.join(prokka_cds.keys())
    # This ensures partial hits in Prokka i.e. Rv####_1, etc is considered
    if locus_tag not in prokka_cds.keys() and locus_tag in prokka_cds_keys:
        found_corresponding_prokka_annotation = True
        same_location = False
        return_feature = None
        partial_annotations = {}
        partial_annotation_keys = {}
        for prokka_locus_tag in prokka_cds.keys():
            if locus_tag in prokka_locus_tag:
                partial_feature = prokka_cds[prokka_locus_tag]
                for feature in partial_feature:
                    feature_key = prokka_locus_tag + ':' + str(feature.location.start) + ':' + \
                                  str(feature.location.end) + ':' + str(feature.location.strand)
                    ann_feature_key = locus_tag + ':' + str(feature.location.start) + ':' + str(feature.location.end) +\
                                      ':' + str(feature.location.strand)
                    partial_annotation_keys[ann_feature_key] = feature_key
                    if feature_key not in partial_feature:
                        partial_annotations[feature_key] = [feature]
                    else:
                        partial_annotations[feature_key].append(feature)
        ratt_feature = locus_tag + ':' + str(ratt_gene_start) + ':' + str(ratt_gene_end) + ':' + str(ratt_gene_strand)
        if ratt_feature in partial_annotation_keys.keys():
            same_location = True
            prokka_partial_annotation_key = partial_annotation_keys[ratt_feature]
            return_feature_dict[prokka_partial_annotation_key] = partial_annotations[prokka_partial_annotation_key]
        else:
            return_feature_dict = partial_annotations
    elif locus_tag in prokka_cds.keys():
        found_corresponding_prokka_annotation = True
        prokka_features_for_gene = prokka_cds[locus_tag]
        same_location = False
        for feature in prokka_features_for_gene:
            feature_start = int(feature.location.start)
            feature_end = int(feature.location.end)
            feature_strand = int(feature.location.strand)
            if feature_strand == ratt_gene_strand and feature_start == ratt_gene_start and feature_end == ratt_gene_end:
                same_location = True
                break
        if same_location:
            return_feature = None
            return_feature_dict = None
        else:
            return_feature = prokka_features_for_gene
    else:
        found_corresponding_prokka_annotation = False
        same_location = False
        return_feature = None
        return_feature_dict = None
    if return_feature is not None:
        for feature in return_feature:
            feature_key = locus_tag + ':' + str(feature.location.start) + ':' + str(feature.location.end) + ':' + \
                          str(feature.location.strand)
            return_feature_dict[feature_key] = feature
    return found_corresponding_prokka_annotation, same_location, return_feature_dict


def get_nonoverlapping_ratt_features(ratt_features_list, prokka_cds_dictionary):
    all_cds = {}
    checked_genes = {}
    overlapping_features = []
    valid_cds = []
    prokka_features_to_add = {}
    genes_in_positive_strand = get_sorted_strand_genes(ratt_features_list, +1)
    genes_in_negative_strand = get_sorted_strand_genes(ratt_features_list, -1)
    for feature in ratt_features_list:
        if feature.type != 'CDS':
            continue
        if len(all_cds.keys()) == 0:
            all_cds[feature.qualifiers['locus_tag'][0]] = [feature]
        elif feature.qualifiers['locus_tag'][0] not in all_cds.keys():
            all_cds[feature.qualifiers['locus_tag'][0]] = [feature]
        # Accounting for duplicate CDSs
        else:
            all_cds[feature.qualifiers['locus_tag'][0]].append(feature)
    pos_strand_keys = genes_in_positive_strand.keys()
    for start in range(1, len(pos_strand_keys)):
        previous_cds_start = pos_strand_keys[start-1]
        current_cds_start = pos_strand_keys[start]
        if previous_cds_start <= current_cds_start <= genes_in_positive_strand[previous_cds_start][0]:
            element_1 = genes_in_positive_strand[previous_cds_start][1] + ':' + str(previous_cds_start) + ':' + \
                        str(genes_in_positive_strand[previous_cds_start][0]) + ':1'
            element_2 = genes_in_positive_strand[current_cds_start][1] + ':' + str(current_cds_start) + ':' + \
                        str(genes_in_positive_strand[current_cds_start][0]) + ':1'
            overlapping_features.append((element_1, element_2))
    neg_strand_keys = genes_in_negative_strand.keys()
    for start in range(1, len(neg_strand_keys)):
        previous_cds_start = neg_strand_keys[start-1]
        current_cds_start = neg_strand_keys[start]
        if previous_cds_start <= current_cds_start <= genes_in_negative_strand[previous_cds_start][0]:
            element_1 = genes_in_negative_strand[previous_cds_start][1] + ':' + str(previous_cds_start) + ':' + \
                        str(genes_in_negative_strand[previous_cds_start][0]) + ':-1'
            element_2 = genes_in_negative_strand[current_cds_start][1] + ':' + str(current_cds_start) + ':' + \
                        str(genes_in_negative_strand[current_cds_start][0]) + ':-1'
            overlapping_features.append((element_1, element_2))
    for overlapping_gene_pair in overlapping_features:
        overlapping_features_elements_1 = overlapping_gene_pair[0].split(':')
        overlapping_features_elements_2 = overlapping_gene_pair[1].split(':')
        gene_1 = overlapping_features_elements_1[0]
        start_gene_1 = int(overlapping_features_elements_1[1])
        end_gene_1 = int(overlapping_features_elements_1[2])
        strand = int(overlapping_features_elements_1[3])  # Strand is the same for both genes
        gene_2 = overlapping_features_elements_2[0]
        start_gene_2 = int(overlapping_features_elements_2[1])
        end_gene_2 = int(overlapping_features_elements_2[2])
        # Adding these genes to checked genes dictionary
        if gene_1 not in checked_genes.keys():
            checked_genes[gene_1] = [':'.join([str(start_gene_1), str(end_gene_1), str(strand)])]
        else:
            checked_genes[gene_1].append(':'.join([str(start_gene_1), str(end_gene_1), str(strand)]))
        if gene_2 not in checked_genes.keys():
            checked_genes[gene_2] = [':'.join([str(start_gene_2), str(end_gene_2), str(strand)])]
        else:
            checked_genes[gene_2].append(':'.join([str(start_gene_2), str(end_gene_2), str(strand)]))
        # Check if gene is possibly merged
        if start_gene_1 == start_gene_2 and end_gene_1 == end_gene_2:
            # If genes are merged, add both the genes in merged position to valid CDS list (possible duplicates in valid
            #  CDS list)
            features_gene_1 = all_cds[gene_1]
            if len(features_gene_1) > 1:
                for feature in features_gene_1:
                    if int(feature.location.start) == start_gene_1 and int(feature.location.end) == end_gene_1 and int(
                            feature.location.strand) == strand:
                        valid_cds.append(feature)
            else:
                valid_cds.append(features_gene_1[0])
            features_gene_2 = all_cds[gene_2]
            if len(features_gene_2) > 1:
                for feature in features_gene_2:
                    if int(feature.location.start) == start_gene_2 and int(feature.location.end) == end_gene_2 and int(
                            feature.location.strand) == strand:
                        valid_cds.append(feature)
                else:
                    valid_cds.append(features_gene_2[0])
            continue
        # Check if neither of the overlapping genes were corrected by RATT
        if gene_1 not in ratt_corrected_genes and gene_2 not in ratt_corrected_genes:
            added_gene_1 = False
            added_gene_2 = False
            # If no correction was performed, these are considered valid overlapping genes
            features_gene_1 = all_cds[gene_1]
            if len(features_gene_1) > 1:
                for feature in features_gene_1:
                    if int(feature.location.start) == start_gene_1 and int(feature.location.end) == end_gene_1 and int(
                            feature.location.strand) == strand:
                        added_gene_1 = True
                        valid_cds.append(feature)
            else:
                added_gene_1 = True
                valid_cds.append(features_gene_1[0])
            features_gene_2 = all_cds[gene_2]
            if len(features_gene_2) > 1:
                for feature in features_gene_2:
                    if int(feature.location.start) == start_gene_2 and int(feature.location.end) == end_gene_2 and int(
                            feature.location.strand) == strand:
                        added_gene_2 = True
                        valid_cds.append(feature)
            else:
                added_gene_2 = True
                valid_cds.append(features_gene_2[0])
            continue
        if gene_1 in ratt_corrected_genes:
            prokka_has_annotation, same_location, corresponding_prokka_feature = \
                check_corresponding_prokka_annotation((gene_1, start_gene_1, end_gene_1, strand), prokka_cds_dictionary)
            if prokka_has_annotation and corresponding_prokka_feature is None:
                # This means that Prokka and Ratt annotated the gene in the same location
                features_gene_1 = all_cds[gene_1]
                if len(features_gene_1) > 1:
                    for feature in features_gene_1:
                        if int(feature.location.start) == start_gene_1 and int(feature.location.end) == end_gene_1 \
                                and int(feature.location.strand) == strand:
                            valid_cds.append(feature)
                else:
                    valid_cds.append(features_gene_1[0])
            elif prokka_has_annotation and corresponding_prokka_feature is not None:
                # This means that Prokka and Ratt annotated the gene in the different locations and we are including
                # annotation from Prokka
                for feature_key in corresponding_prokka_feature.keys():
                    if feature_key not in prokka_features_to_add.keys():
                        prokka_features_to_add[feature_key] = [corresponding_prokka_feature[feature_key]]
                    else:
                        for feature in corresponding_prokka_feature[feature_key]:
                            prokka_features_to_add[feature_key].append(feature)
            elif not prokka_has_annotation and corresponding_prokka_feature is None:
                # This means that Prokka annotation is missing for corresponding gene. Ratt annotation is taken and
                # feature printed to screen
                features_gene_1 = all_cds[gene_1]
                if len(features_gene_1) > 1:
                    for feature in features_gene_1:
                        if int(feature.location.start) == start_gene_1 and int(feature.location.end) == end_gene_1 \
                                and int(feature.location.strand) == strand:
                            valid_cds.append(feature)
                else:
                    valid_cds.append(features_gene_1[0])
                print('The gene ' + gene_1 + ' overlaps with ' + gene_2 + ' and has been corrected by RATT but this '
                                                                          'gene is not annotated by Prokka so RATT '
                                                                          'annotation is kept')
                print(overlapping_gene_pair)
            ### Can Be Removed After Debugging ###
            else:
                print('Corner Case!!!!!!!!!!!')
                print(overlapping_gene_pair)
                print(prokka_has_annotation)
                print(corresponding_prokka_feature)
            #######################################
        if gene_2 in ratt_corrected_genes:
            prokka_has_annotation, same_location, corresponding_prokka_feature = \
                check_corresponding_prokka_annotation((gene_2, start_gene_2, end_gene_2, strand), prokka_cds_dictionary)
            if prokka_has_annotation and corresponding_prokka_feature is None:
                # This means that Prokka and Ratt annotated the gene in the same location
                features_gene_2 = all_cds[gene_2]
                if len(features_gene_2) > 1:
                    for feature in features_gene_2:
                        if int(feature.location.start) == start_gene_2 and int(feature.location.end) == end_gene_2 \
                                and int(feature.location.strand) == strand:
                            valid_cds.append(feature)
                else:
                    valid_cds.append(features_gene_2[0])
            elif prokka_has_annotation and corresponding_prokka_feature is not None:
                # This means that Prokka and Ratt annotated the gene in the different locations and we are including
                # annotation from Prokka
                for feature_key in corresponding_prokka_feature.keys():
                    if feature_key not in prokka_features_to_add.keys():
                        prokka_features_to_add[feature_key] = [corresponding_prokka_feature[feature_key]]
                    else:
                        for feature in corresponding_prokka_feature[feature_key]:
                            prokka_features_to_add[feature_key].append(feature)
            elif not prokka_has_annotation and corresponding_prokka_feature is None:
                # This means that Prokka annotation is missing for corresponding gene. Ratt annotation is taken and
                # feature printed to screen
                features_gene_2 = all_cds[gene_2]
                if len(features_gene_2) > 1:
                    for feature in features_gene_2:
                        if int(feature.location.start) == start_gene_2 and int(feature.location.end) == end_gene_2 \
                                and int(feature.location.strand) == strand:
                            valid_cds.append(feature)
                else:
                    valid_cds.append(features_gene_2[0])
                print('The gene ' + gene_2 + ' overlaps with ' + gene_1 + ' and has been corrected by RATT but this '
                                                                          'gene is not annotated by Prokka so RATT '
                                                                          'annotation is kept')
                print(overlapping_gene_pair)
            ### Can Be Removed After Debugging ###
            else:
                print('Corner Case!!!!!!!!!!!')
                print(overlapping_gene_pair)
                print(prokka_has_annotation)
                print(corresponding_prokka_feature)
            #######################################
    # Add all 'non-overlapping' gene instances to valid RATT cds list
    for gene in all_cds.keys():
        gene_features = all_cds[gene]
        if gene not in checked_genes.keys():
            for feature in gene_features:
                valid_cds.append(feature)
        else:
            checked_positions = checked_genes[gene]
            for feature in gene_features:
                feature_key = ':'.join([str(feature.location.start), str(feature.location.end),
                                        str(feature.location.strand)])
                if feature_key not in checked_positions:
                    valid_cds.append(feature)
    return valid_cds, prokka_features_to_add


def genes_locus_tag_parser(file_path=None):
    """
    :param file_path: Filepath for genes.tsv where 1st column is Rv locus tag and second column is gene name if present.
    Default file used is $GROUPHOME/resources/mtb-reconstruction/genes.tsv
    :return: Dictionary where keys are gene names and values are corresponding Rv locus tags.
    """
    genes_dict = {}
    if file_path is None:
        file_path = os.environ['GROUPHOME'] + '/resources/mtb-reconstruction/genes.tsv'
    genes_raw = open(file_path, 'r').readlines()
    for line in genes_raw:
        if line[0] != 'R':
            continue
        line_elements = line.strip().split()
        if line_elements[1] == line_elements[0]:
            continue
        else:
            genes_dict[line_elements[1]] = line_elements[0]
    return genes_dict


def get_prokka_cds(all_annotations_prokka, only_rv=False):
    prokka_cds_dict = {}
    global gene_name_to_rv
    gene_name_to_rv = genes_locus_tag_parser()
    if only_rv:
        for feature in all_annotations_prokka:
            if feature.type != 'CDS' or 'gene' not in feature.qualifiers.keys():
                continue
            else:
                dict_key = feature.qualifiers['gene'][0]
                if dict_key[:2] != 'Rv':
                    if dict_key in gene_name_to_rv.keys():
                        locus_tag = gene_name_to_rv[dict_key]
                    else:
                        locus_tag = dict_key
                else:
                    locus_tag = dict_key
                if locus_tag not in prokka_cds_dict.keys():
                    prokka_cds_dict[locus_tag] = [feature]
                else:
                    prokka_cds_dict[locus_tag].append(feature)
    else:
        for feature in all_annotations_prokka:
            if feature.type != 'CDS':
                continue
            if 'gene' in feature.qualifiers.keys():
                dict_key = feature.qualifiers['gene'][0]
            else:
                dict_key = feature.qualifiers['locus_tag'][0]
            if dict_key[:2] != 'Rv':
                if dict_key in gene_name_to_rv.keys():
                    locus_tag = gene_name_to_rv[dict_key]
                else:
                    locus_tag = dict_key
            else:
                locus_tag = dict_key
            if locus_tag not in prokka_cds_dict.keys():
                prokka_cds_dict[locus_tag] = [feature]
            else:
                prokka_cds_dict[locus_tag].append(feature)
    return prokka_cds_dict


def blast_feature_sequence_to_rv(query, locus_tag, to_print=False):
    #add_annotation = False
    if locus_tag in h37rv_sequences.keys():
        subject_fp = '/tmp/' + locus_tag + '.fasta'
        blast_to_rv = NcbiblastpCommandline(subject=subject_fp, outfmt='"7 qseqid qlen sseqid slen qlen length'
                                                                               ' pident qcovs"')
        stdout, stderr = blast_to_rv(stdin=query)
        if to_print:
            print(stdout)
        hit, all_hits, note = identify_top_hits(stdout)
        if hit == locus_tag:
            add_annotation = True
        else:
            add_annotation = False
    else:
        add_annotation = True
#        print('No sequence found for locus tag in fasta file: ' + locus_tag)
    return add_annotation


def isolate_valid_ratt_annotations(feature_list, prokka_annotations):
    # This function takes as input a list of features and checks if the length of the CDSs are divisible by 3. If not,
    # it outputs the features to stdout and removes them from the valid_ratt_annotations
    global valid_amino_acids
    valid_amino_acids = ['M', 'V', 'L']
    unbroken_cds = []
    genes_from_prokka = {}
    num_joins = 0
    non_cds_features = []
    prokka_cds_to_add = []
    prokka_added = {}
    broken_cds = []
    for feature in feature_list:
        if feature.type == 'CDS' and (len(feature.location) % 3) != 0:
            print('Nucleotide sequence is not divisible by 3')
            print(feature)
        # Identify features with 'joins'
        elif feature.type == 'CDS' and 'Bio.SeqFeature.CompoundLocation' in str(type(feature.location)):
            locus_tag = feature.qualifiers['locus_tag'][0]
            gene_name = locus_tag
            if locus_tag not in broken_cds:
                broken_cds.append(locus_tag)
            if 'gene' in feature.qualifiers.keys():
                gene_name = feature.qualifiers['gene'][0]
            if locus_tag not in genes_from_prokka.keys():
                genes_from_prokka[locus_tag] = gene_name
            num_joins += 1
        elif feature.type == 'CDS':
            unbroken_cds.append(feature)
        else:
            non_cds_features.append(feature)
    prokka_cds_dict = get_prokka_cds(prokka_annotations, only_rv=True)
    valid_cds, to_add_from_prokka = get_nonoverlapping_ratt_features(unbroken_cds, prokka_cds_dict)
    valid_features = []
    print("Valid CDSs before checking coverage: " + str(len(valid_cds)))
    for cds_feature in valid_cds:
        feature_sequence = translate(cds_feature.extract(record_sequence), to_stop=True)
        cds_locus_tag = cds_feature.qualifiers['locus_tag'][0]
        if len(feature_sequence) == 0 or feature_sequence[0] not in valid_amino_acids:
            continue
        else:
            add_sequence = blast_feature_sequence_to_rv(str(feature_sequence), cds_locus_tag)
            if add_sequence:
                cds_feature.qualifiers['translation'] = [str(feature_sequence)]
                valid_features.append(cds_feature)
            else:
                continue

    print("Valid CDSs after checking coverage: " + str(len(valid_features)))
    for non_cds in non_cds_features:
        valid_features.append(non_cds)
    for gene_key in to_add_from_prokka.keys():
        prokka_cds_to_add.append(to_add_from_prokka[gene_key][0])
        locus_tag = gene_key.split(':')[0]
        prokka_added[locus_tag] = [gene_key]
    for gene in genes_from_prokka.keys():
        prokka_gene_features = {}
        gene_name = genes_from_prokka[gene]
        if gene_name in prokka_cds_dict.keys():
            prokka_gene_features[gene] = prokka_cds_dict[gene_name]
        elif gene in prokka_cds_dict.keys():
            prokka_gene_features[gene] = prokka_cds_dict[gene]
        else:
            continue
        for feature in prokka_gene_features[gene]:
            feature_key = ':'.join([gene, str(feature.location.start), str(feature.location.end),
                                    str(feature.location.strand)])
            if gene in prokka_added.keys() and feature_key in prokka_added[gene]:
                continue
            elif gene in prokka_added.keys() and feature_key not in prokka_added[gene]:
                prokka_cds_to_add.append(feature)
                prokka_added[gene].append(feature_key)
            else:
                prokka_cds_to_add.append(feature)
                prokka_added[gene] = [feature_key]
    return valid_features, prokka_cds_to_add


def remove_duplicate_annotations(ratt_features, prokka_features_dictionary):
    # This function prunes and selects the features that are relevant in Prokka and
    # discards features in Prokka that are annotated by RATT by taking into account the
    # gene name and the position of the features
    prokka_features_not_in_ratt = prokka_features_dictionary.copy()
    ratt_overlapping_genes = {}
    for ratt_feature in ratt_features:
        if ratt_feature.type == 'gene' or ratt_feature.type == 'mRNA':
            continue
        ratt_start = int(ratt_feature.location.start)
        ratt_end = int(ratt_feature.location.end)
        ratt_strand = ratt_feature.location.strand
        prokka_duplicate_removed = False
        for prokka_feature_position in prokka_features_dictionary.keys():
            prokka_feature = prokka_features_dictionary[prokka_feature_position]
            prokka_start = prokka_feature_position[0]
            prokka_end = prokka_feature_position[1]
            prokka_strand = prokka_feature_position[2]
            if prokka_start > ratt_start and prokka_end > ratt_start:
                break
            elif (abs(ratt_start-prokka_start) <= 5) and (abs(ratt_end-prokka_end) <= 5):
                if ('gene' in ratt_feature.qualifiers.keys() and 'gene' in prokka_feature.qualifiers.keys()) and \
                        (ratt_feature.qualifiers['gene'] == prokka_feature.qualifiers['gene']):
                    prokka_features_not_in_ratt.pop((prokka_start,prokka_end,prokka_strand), None)
                    prokka_duplicate_removed = True
                elif len(ratt_feature.location) == len(prokka_feature.location):
                    prokka_features_not_in_ratt.pop((prokka_start,prokka_end,prokka_strand), None)
                    prokka_duplicate_removed = True
        if not prokka_duplicate_removed and ratt_feature.type != 'mRNA':
            ratt_overlapping_genes[(ratt_start, ratt_end, ratt_strand)] = ratt_feature
    return prokka_features_not_in_ratt, ratt_overlapping_genes


#################################################################################
# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
# This function was modified by Deepika Gunasekaran
def get_interregions(embl_record,intergene_length=1):
    seq_record = embl_record
    cds_list_plus = []
    cds_list_minus = []
    intergenic_records = []
    intergenic_positions = []
    pre_intergene = {}
    post_intergene = {}
    # Loop over the genome file, get the CDS features on each of the strands
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        mystart = feature.location.start.position
        myend = feature.location.end.position
        if feature.strand == -1:
            cds_list_minus.append((mystart,myend,-1))
            pre_intergene[(myend,-1)] = feature
            post_intergene[(mystart,-1)] = feature
        elif feature.strand == 1:
            cds_list_plus.append((mystart,myend,1))
            pre_intergene[(myend,1)] = feature
            post_intergene[(mystart,1)] = feature
        else:
            sys.stderr.write("No strand indicated %d-%d. Assuming +\n" %
                                  (mystart, myend))
            cds_list_plus.append((mystart,myend,1))
            pre_intergene[(myend,1)] = feature
            post_intergene[(mystart,1)] = feature
    cds_list_plus = sorted(cds_list_plus)
    cds_list_minus = sorted(cds_list_minus)
    for i,pospair in enumerate(cds_list_plus[1:]):
        # Compare current start position to previous end position
        last_end = cds_list_plus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end:this_start]
            strand_string = "+"
            intergenic_records.append(
                  SeqRecord(intergene_seq,id="%s-ign-%d" % (seq_record.name,i),
                  description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                        this_start,strand_string)))
            intergenic_positions.append((last_end+1,this_start,strand_string))
    for i,pospair in enumerate(cds_list_minus[1:]):
        last_end = cds_list_minus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end:this_start]
            strand_string = "-"
            intergenic_records.append(
                  SeqRecord(intergene_seq,id="%s-ign-%d" % (seq_record.name,i),
                  description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                        this_start,strand_string)))
            intergenic_positions.append((last_end+1,this_start,strand_string))
    return intergenic_records, intergenic_positions, pre_intergene, post_intergene
#####################################################################################


def check_inclusion_criteria(annotation_mapping_dict, embl_file, ratt_annotation, prokka_annotation):
    ratt_feature_start = int(ratt_annotation.location.start)
    ratt_feature_end = int(ratt_annotation.location.end)
    prokka_feature_start = int(prokka_annotation.location.start)
    prokka_feature_end = int(prokka_annotation.location.end)
    included = False
    # Check if feature types are the same. If not add feature to EMBL record
    if ratt_annotation.type != prokka_annotation.type:
        if prokka_annotation.type == 'CDS':
            mod_prokka_annotation = validate_prokka_feature_annotation(prokka_annotation, prokka_noref_dictionary)
        else:
            mod_prokka_annotation = prokka_annotation
        embl_file.features.append(mod_prokka_annotation)
        included = True
    # Check if gene names match and if they don't or if gene names are missing, keep both
    elif 'gene' in ratt_annotation.qualifiers.keys() and 'gene' in prokka_annotation.qualifiers.keys():
        if ratt_annotation.qualifiers['gene'] != prokka_annotation.qualifiers['gene']:
            mod_prokka_annotation = validate_prokka_feature_annotation(prokka_annotation, prokka_noref_dictionary)
            embl_file.features.append(mod_prokka_annotation)
            included = True
        # If gene names are the same and the lengths of the genes are comparable between RATT and Prokka annotation
        # (difference in length of less than/equal to 10 bps), the RATT annotation is prefered
        elif ratt_annotation.qualifiers['gene'] == prokka_annotation.qualifiers['gene'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            mod_prokka_annotation = validate_prokka_feature_annotation(prokka_annotation, prokka_noref_dictionary)
            embl_file.features.append(mod_prokka_annotation)
            included = True
        # If gene tag is missing and the product is not a hypothetical protein, check to see if the products are the
            # same between RATT and Prokka and if they are and if the lengths of the protein coding genes are comparable
            # between RATT and Prokka (difference in length is less than/equal to 10 bps), the RATT annotation is
            # prefered
    elif ('gene' not in ratt_annotation.qualifiers.keys() or 'gene' not in prokka_annotation.qualifiers.keys()) and \
            ('product' in ratt_annotation.qualifiers.keys() and 'product' in prokka_annotation.qualifiers.keys()):
        if ratt_annotation.qualifiers['product'] == 'hypothetical protein' or \
                        prokka_annotation.qualifiers['product'] == 'hypothetical protein':
            mod_prokka_annotation = validate_prokka_feature_annotation(prokka_annotation, prokka_noref_dictionary)
            embl_file.features.append(mod_prokka_annotation)
            included = True
        elif ratt_annotation.qualifiers['product'] == prokka_annotation.qualifiers['product'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            mod_prokka_annotation = validate_prokka_feature_annotation(prokka_annotation, prokka_noref_dictionary)
            embl_file.features.append(mod_prokka_annotation)
            included = True
    return embl_file, included


def get_top_hit(all_hits_dict):
    """
    :param all_hits_dict: This dictionary consists of all Blast hits with corresponding unannotated locus_tag that have
    passed the defined thresholds for amino acid identity and coverage
    :return: returns only the elements of top hits with the locus tag. The key is the 'L(2)_' locus tag and the value is
     corresponding top hit in H37Rv. If multiple top hits are present with same identity and coverage values, the value
     is all the top Rv hits separated by ':'
    """
    top_hit_identity = 0
    for gene in all_hits_dict.keys():
        if all_hits_dict[gene] > top_hit_identity:
            top_hit_identity = all_hits_dict[gene]
            top_hit = gene
    return top_hit


def identify_top_hits(blast_output_file, identity=95, coverage=95, mtb=False):
    """
    :param blast_output_file: Output file from Blastp
    :return: Dictionary of top hits for each unannotated locus that passes the identity and coverage threshold
    """
    if mtb:
        all_hits_dict = {}
        blast_output_raw = blast_output_file.split('\n')
        mtb_hit = False
        for line in blast_output_raw:
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            line_elements = line.split('\t')
            iden = float(line_elements[5])
            qcov = (float(line_elements[4]) / float(line_elements[1])) * 100.0
            scov = (float(line_elements[4]) / float(line_elements[3])) * 100.0
            if iden >= identity and qcov >= coverage and scov >= coverage:
                mtb_hit = True
                corresponding_mtb_hit = line_elements[2]
                all_hits_dict[corresponding_mtb_hit] = iden
        if mtb_hit:
            top_hit = get_top_hit(all_hits_dict)
        else:
            top_hit = None
        return top_hit, all_hits_dict
    else:
        all_hits_dict = {}
        notes = {}
        blast_output_raw = blast_output_file.split('\n')
        rv_hit = False
        for line in blast_output_raw:
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            line_elements = line.split('\t')
            iden = float(line_elements[5])
            qcov = (float(line_elements[4]) / float(line_elements[1])) * 100.0
            scov = (float(line_elements[4]) / float(line_elements[3])) * 100.0
            if iden >= identity and qcov >= coverage and scov >= coverage:
                rv_hit = True
                corresponding_rv_hit = line_elements[2].split('|')[-1]
                all_hits_dict[corresponding_rv_hit] = iden
            elif 85 <= iden < identity and qcov >= coverage and scov >= coverage:
                corresponding_rv_hit = line_elements[2].split('|')[-1]
                notes[corresponding_rv_hit] = iden
            else:
                continue
        if rv_hit:
            top_hit = get_top_hit(all_hits_dict)
        elif len(all_hits_dict.keys()) != 0:
            top_hit = 'MTB:' + get_top_hit(all_hits_dict)
        else:
            top_hit = None
        return top_hit, all_hits_dict, notes


def get_essential_genes():
    gene_essentiality_fp = os.environ['GROUPHOME'] +'/resources/gene_essentiality_dejesus_2017.tsv'
    gene_essentiality_raw = open(gene_essentiality_fp, 'r').readlines()
    essential_genes = []
    for line in gene_essentiality_raw:
        line_elements = line.strip().split('\t')
        if line_elements[1] == 'ES':
            essential_genes.append(line_elements[0])
    return essential_genes


def get_lengths_of_genes(genes_list, fasta_fp, fasta_conversion=False):
    fasta_dict = {}
    fasta_records = SeqIO.parse(fasta_fp, 'fasta')
    for record in fasta_records:
        header_raw = record.description
        header_elements = header_raw.strip().split(' [locus_tag=')[1]
        rv_id = header_elements.split(']')[0]
        sequence = str(record.seq)
        if not fasta_conversion:
            fasta_dict[rv_id] = len(sequence)
        else:
            fasta_dict[rv_id] = len(sequence) * 3
    return fasta_dict


def get_rv_sequences_from_fasta(h37rv_fasta_fp):
    rv_sequence_dict = {}
    fasta_records = SeqIO.parse(h37rv_fasta_fp, 'fasta')
    for record in fasta_records:
        header_raw = record.description
        header_elements = header_raw.strip().split(' [locus_tag=')[1]
        rv = header_elements.split(']')[0]
        rv_seq = str(record.seq)
        fp = '/tmp/' + rv + '.fasta'
        header = '>' + rv + '\n'
        seq = rv_seq
        rv_sequence_dict[rv] = seq
        with open(fp, 'w') as fasta_file:
            fasta_file.write(header)
            fasta_file.write(seq)
    return rv_sequence_dict


def validate_prokka_feature_annotation(feature, prokka_noref):
    mod_feature = feature
    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
    if feature.type != 'CDS':
        return mod_feature
    if 'gene' not in feature.qualifiers.keys():
        mod_feature = feature
    else:
        feature_seq = str(feature.qualifiers['translation'][0])
        if feature.qualifiers['gene'][0][:2] == 'Rv':
            locus_tag = feature.qualifiers['gene'][0]
#            print('Blasting ' + locus_tag + 'to H37Rv')
            prokka_blast_list.append(locus_tag)
            retain_rv_annotation = blast_feature_sequence_to_rv(feature_seq, locus_tag)
            if retain_rv_annotation:
                mod_feature = feature
            else:
                mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
        elif feature.qualifiers['gene'][0] in gene_name_to_rv.keys():
            locus_tag = gene_name_to_rv[feature.qualifiers['gene'][0]]
#            print('Blasting ' + locus_tag + 'to H37Rv')
            prokka_blast_list.append(locus_tag)
            retain_rv_annotation = blast_feature_sequence_to_rv(feature_seq, locus_tag)
            if retain_rv_annotation:
                mod_feature = feature
            else:
                mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
        else:
            mod_feature = prokka_noref[loc_key]
    return mod_feature


# Required inputs:
# 1. Isolate ID
# Optional inputs:
# 1. Output file/path for Genbank output
# 2. Output file/path for log file


def main():
    input_ratt_files = []
    # Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file> -m <output_merged_genes>
    # from the commandline. The log file output stats about the features that are added to the RATT annotation. The
    # default locations for the -g and -l options are 'isolate_id'/annomerge/'isolate_id'.gbf and
    # 'isolate_id'/annomerge/'isolate_id'.log
    parser = argparse.ArgumentParser(description='Merging annotation from RATT and Prokka',
                                     epilog='Isolate ID must be specified')
    parser.add_argument('-i', '--isolate', help='Isolate ID')
    parser.add_argument('-o', '--output', help='Output file in Genbank format', default='annomerge.gbf')
    parser.add_argument('-l', '--log_file', help='Log file with information on features added from prokka',
                        default='annomerge.log')
    parser.add_argument('-m', '--merged_genes', help='Merged genes file in genbank format',
                        default='merged_genes.gbf')
    parser.add_argument('-illumina', '--illumina', help='Set this flag to true if isolate has no circularized contigs '
                                                        'or if dnaA need not be verified',
                        default=False)
    args = parser.parse_args()

    if args.isolate is None:
        parser.print_help()
        sys.exit('Isolate ID must be specified')
    file_path = os.environ['GROUPHOME'] + '/data/depot/annotation/' + args.isolate + '/'
    ratt_file_path = file_path + 'ratt'
    ratt_correction_files = []
    try:
        ratt_embl_files = [embl_file for embl_file in os.listdir(ratt_file_path) if embl_file.endswith('.final.embl')]
        correction_files = [cf for cf in os.listdir(ratt_file_path) if cf.endswith('.Report.txt')]
        for embl_file in ratt_embl_files:
            embl_file_path = ratt_file_path + '/' + embl_file
            input_ratt_files.append(embl_file_path)
        for corr_file in correction_files:
            corr_file_path = ratt_file_path + '/' + corr_file
            ratt_correction_files.append(corr_file_path)
    except OSError:
        sys.exit('Expecting RATT annotation files but found none')
    try:
        input_prokka_genbank = file_path + 'prokka/' + args.isolate + '.gbf'
    except OSError:
        sys.exit('Expecting Prokka annotation file but found none')
    output_merged_genes = args.merged_genes
    output_genbank = args.output
    output_log = args.log_file
    global output_file
    output_file = open(output_log, 'w+')
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    prokka_record_fp = file_path + 'prokka-noreference/' + args.isolate + '.gbf'
    prokka_record_noref = list(SeqIO.parse(prokka_record_fp, 'genbank'))
    annomerge_records = []
    h37rv_protein_fasta = os.environ['GROUPHOME'] + '/resources/H37Rv-CDS-fullheader-NC_000962.3.fasta'
    global h37rv_sequences
    h37rv_sequences = get_rv_sequences_from_fasta(h37rv_protein_fasta)
    global h37rv_protein_lengths
    global genes_in_rv
    genes_in_rv = []
    genes_in_rv_file = os.environ['GROUPHOME'] + '/resources/H37rv_proteincoding_genes.txt'
    genes_in_rv_raw = open(genes_in_rv_file, 'r').readlines()
    for line in genes_in_rv_raw:
        locus_tag = line.split()[0]
        genes_in_rv.append(locus_tag)
    h37rv_protein_lengths = get_lengths_of_genes(genes_in_rv, h37rv_protein_fasta)
#    print(input_ratt_files)
#    print(len(prokka_records))
    for i in range(0, len(input_ratt_files)):
        ratt_contig_record = SeqIO.read(input_ratt_files[i], 'embl')
        prokka_noref_rec = prokka_record_noref[i]
        global prokka_noref_dictionary
        prokka_noref_dictionary = generate_feature_dictionary(prokka_noref_rec.features)
        global record_sequence
        record_sequence = ratt_contig_record.seq
        prokka_contig_record = prokka_records[i]
        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        prokka_rv_features, prokka_non_rv_features = get_prokka_features_with_gene_name(prokka_contig_features)
        if i == 0 and not args.illumina:
            check_for_dnaA(ratt_contig_features)
        if args.illumina:
            print('DnaA might not be the first element. Circularization is not checked')
        if len(ratt_correction_files) == 1:
            error_correction_fp = ratt_correction_files[0]
        else:
            error_correction_fp = ratt_correction_files[i]
        global ratt_corrected_genes
        ratt_corrected_genes = get_ratt_corrected_genes(error_correction_fp)
        ratt_contig_features, prokka_features_to_add = isolate_valid_ratt_annotations(ratt_contig_features,
                                                                                      prokka_contig_features)
        ratt_contig_features = remove_duplicate_cds(ratt_contig_features)
        cds_from_ratt = 0
        for feature in ratt_contig_features:
            if feature.type == 'CDS':
                cds_from_ratt += 1
        ratt_contig_features = get_ordered_features(ratt_contig_features)
        merged_genes, check_prokka = identify_merged_genes(ratt_contig_features)
        if check_prokka:
            merged_features, corner_cases, corner_cases_explicit, ratt_contig_features_mod, \
            prokka_contig_features_mod = get_annotation_for_merged_genes(merged_genes, prokka_contig_features,
                                                                         ratt_contig_features)
            ratt_contig_features = ratt_contig_features_mod
            prokka_contig_features = prokka_contig_features_mod
            merged_features_record = prokka_contig_record[:]
            merged_features_record.features = merged_features
            if corner_cases:
                print('MERGED GENES: Corner cases')
                for strand in corner_cases_explicit.keys():
                    if len(corner_cases_explicit[strand]) > 0:
                        print(args.isolate + ' '.join(corner_cases_explicit[strand]))
        else:
            merged_features = []
        if len(merged_features) > 0 and i == 0:
            SeqIO.write(merged_features_record, output_merged_genes, 'genbank')
        ratt_cds_count = 0
        for feature in ratt_contig_features:
            if feature.type == 'CDS':
                ratt_cds_count += 1
#        output_file.write('Number of valid CDSs in RATT: ' + str(ratt_cds_count) + '\n')
#        output_file.write('Number of merged genes in contig ' + str(i + 1) + ' = ' + str(len(merged_features)) + '\n')
        if len(ratt_contig_features) == 0:
            print("NO RATT ANNOTATION FOR CONTIG " + str(i + 1))
            feature_additions = {}
            feature_lengths = {}
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
            # print("Annomerge record length")
            # print(len(annomerge_records))
            for prokka_feature in prokka_contig_record.features:
                if prokka_feature.type not in feature_additions.keys():
                    feature_additions[prokka_feature.type] = 1
                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                else:
                    feature_additions[prokka_feature.type] += 1
                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            output_file.write('Contig Number: ' + str(i + 1) + '\n')
            for f in feature_lengths.keys():
                output_file.write(str('Feature: ' + f + '\n'))
                output_file.write(str('Number of ' + f + ': ' + str(feature_additions[f]) + '\n'))
                output_file.write(str('Min length: ' + str(min(feature_lengths[f])) + '\n'))
                output_file.write(str('Max length: ' + str(max(feature_lengths[f])) + '\n'))
                output_file.write(str('Median length: ' + str(median(feature_lengths[f])) + '\n'))
            # annomerge_contig_record = prokka_contig_record
            continue
        elif len(prokka_contig_features) == 0:
            print("NO PROKKA ANNOTATION FOR CONTIG " + str(i + 1))
            # annomerge_contig_record = ratt_contig_record
            prokka_contig_features = ratt_contig_features
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            output_file.write('Contig Number: ' + str(i + 1) + '\n')
            output_file.write('No Annotation to add from Prokka')
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the
            # features from RATT
            annomerge_contig_record = prokka_contig_record[:]
            annomerge_contig_record.features = ratt_contig_features
            num_feat = 0
            for test_f in annomerge_contig_record.features:
                if test_f.type == 'CDS':
                    num_feat += 1
#            print('Number of features before adding Prokka')
#            print(num_feat)
            annomerge_contig_record.annotations['comment'] = 'Merged reference based annotation from RATT and ab ' \
                                                             'initio annotation from Prokka'
            ###########################################################
            ####### Creating a dictionary with feature location #######
            ############### as key and index as value #################
            ###########################################################
            ratt_annotation_mapping = {}  # Used for resolving annotations of overlapping features between RATT and
            # Prokka
            for index, feature in enumerate(ratt_contig_record.features):
                start = int(feature.location.start)
                end = int(feature.location.end)
                ratt_annotation_mapping[(start, end)] = index
            prokka_annotation_mapping = {}
            ratt_contig_record_mod = ratt_contig_record[:]
            ratt_contig_record_mod.features = ratt_contig_features
#            print(len(ratt_contig_record_mod.features))
            added_from_ratt = []
            global prokka_blast_list
            prokka_blast_list = []
            for feat in ratt_contig_record_mod.features:
                if feat.type == 'CDS':
                    added_from_ratt.append(feat.qualifiers['locus_tag'][0])
            prokka_features_dict = generate_feature_dictionary(prokka_contig_features)
            prokka_features_not_in_ratt, ratt_overlapping_genes = \
                remove_duplicate_annotations(ratt_contig_features, prokka_features_dict)
            intergenic_ratt, intergenic_positions, ratt_pre_intergene, ratt_post_intergene = \
                get_interregions(ratt_contig_record_mod, intergene_length=1)
            sorted_intergenic_positions = sorted(intergenic_positions)
            feature_additions = {}
            feature_lengths = {}
            for j in sorted_intergenic_positions:
                # Variable definitions
                ratt_unannotated_region_start = j[0]
                ratt_unannotated_region_end = j[1]
                ratt_strand = j[2]
                for feature_position in prokka_features_not_in_ratt.keys():
                    prokka_feature = prokka_features_not_in_ratt[feature_position]
                    prokka_strand = feature_position[2]
                    prokka_feature_start = feature_position[0]
                    prokka_feature_end = feature_position[1]
                    ratt_unannotated_region_range = range(ratt_unannotated_region_start,
                                                          ratt_unannotated_region_end + 1)
                    if (prokka_strand == -1 and ratt_strand == '-') or (prokka_strand == 1 and ratt_strand == '+'):
                        # If Prokka feature is location before the start of the intergenic region, pop the key from the
                        # dictionary and continue loop
                        if (prokka_feature_start < ratt_unannotated_region_start) and \
                                (prokka_feature_end < ratt_unannotated_region_start):
                            prokka_features_not_in_ratt.pop((prokka_feature_start, prokka_feature_end, prokka_strand),
                                                            None)
                            continue
                        # Else if the prokka feature location is after the end of the intergenic region, break out of
                        # the inner loop
                        elif (prokka_feature_start > ratt_unannotated_region_end) and \
                                (prokka_feature_end > ratt_unannotated_region_end):
                            break
                        # If the PROKKA feature is contained in the RATT feature
                        elif prokka_feature_start in ratt_unannotated_region_range and \
                                prokka_feature_end in ratt_unannotated_region_range:
                            # The if-else condition below is to keep track of the features added from Prokka for the
                            # log file
                            if prokka_feature.type not in feature_additions.keys():
                                feature_additions[prokka_feature.type] = 1
                                feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                            else:
                                feature_additions[prokka_feature.type] += 1
                                feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
                            mod_prokka_feature = validate_prokka_feature_annotation(prokka_feature,
                                                                                    prokka_noref_dictionary)
                            annomerge_contig_record.features.append(mod_prokka_feature)
                        # If the Prokka feature overlaps with two RATT features
                        elif prokka_feature_start < ratt_unannotated_region_start and \
                                prokka_feature_end > ratt_unannotated_region_end:
                            if prokka_feature.type == 'source':  # This is to exclude the source feature as it is
                                # accounted for by RATT
                                break
                            #print('Overlaps with 2 RATT genes')
                            #print(prokka_feature)
                            mod_prokka_feature = validate_prokka_feature_annotation(prokka_feature,
                                                                                    prokka_noref_dictionary)
                            annomerge_contig_record.features.append(mod_prokka_feature)
                            # The if-else condition below is to keep track of the features added from Prokka for the
                            # log file
                            if prokka_feature.type not in feature_additions.keys():
                                feature_additions[prokka_feature.type] = 1
                                feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                            else:
                                feature_additions[prokka_feature.type] += 1
                                feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
                        # If the Prokka feature overlaps with one RATT feature
                        else:
                            if (prokka_feature_start < ratt_unannotated_region_start) and \
                                    (prokka_feature_end in ratt_unannotated_region_range):
                                ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                               prokka_strand)]
                            elif (prokka_feature_start in ratt_unannotated_region_range) and \
                                    prokka_feature_end > ratt_unannotated_region_end:
                                ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                                prokka_strand)]
                            #print('Overlaps with 1 ratt')
                            #print(ratt_overlapping_feature)
                            #print(prokka_feature)
                            annomerge_contig_record, included = check_inclusion_criteria(ratt_annotation_mapping,
                                                                                         annomerge_contig_record,
                                                                                         ratt_overlapping_feature,
                                                                                         prokka_feature)
                            if included:  # To check if Prokka feature passed the inclusion criteria and was integrated
                                # into the EMBL file
                                # The if-else condition below is to keep track of the features added from Prokka for the
                                # log file
                                if prokka_feature.type not in feature_additions.keys():
                                    feature_additions[prokka_feature.type] = 1
                                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                                else:
                                    feature_additions[prokka_feature.type] += 1
                                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            if len(merged_features) > 0:
                for feature_1 in merged_features:
                    annomerge_contig_record.features.append(feature_1)

            unflattened_features = prokka_features_to_add
            flattened_features = []
            for feature_type in unflattened_features:
                if 'Bio.SeqFeature.SeqFeature' in str(type(feature_type)):
                    flattened_features.append(feature_type)
                elif 'list' in str(type(feature_type)) and len(feature_type) > 0:
                    for sub_feature in feature_type:
                        if 'Bio.SeqFeature.SeqFeature' not in str(type(sub_feature)):
                            continue
                        else:
                            flattened_features.append(sub_feature)
                else:
                    continue
#            print('Prokka_features_to_add')
#            print(len(flattened_features))
            for feature in flattened_features:
                mod_feature = validate_prokka_feature_annotation(feature, prokka_noref_dictionary)
                annomerge_contig_record.features.append(mod_feature)
            annomerge_records.append(annomerge_contig_record)

    ###############################################################################################
    ######## Post-processing of genbank file to remove duplicates and rename locus_tag for ########
    ################################### Prokka annotations ########################################
    ###############################################################################################
    prokka_record_fp = file_path + 'prokka-noreference/' + args.isolate + '.gbf'
    prokka_record_noref = list(SeqIO.parse(prokka_record_fp, 'genbank'))
    counter = 0
    annomerge_records_post_processed = []
    for rec_num in range(0, len(annomerge_records)):
        prokka_rec = annomerge_records[rec_num]
        prokka_noref_rec = prokka_record_noref[rec_num]
        prokka_noref_dict = generate_feature_dictionary(prokka_noref_rec.features)
        raw_features_unflattened = prokka_rec.features[:]
        raw_features = []
        for f_type in raw_features_unflattened:
            if 'Bio.SeqFeature.SeqFeature' in str(type(f_type)):
                raw_features.append(f_type)
            elif 'list' in str(type(f_type)) and len(f_type) > 0:
                for sub_feature in f_type:
                    if 'Bio.SeqFeature.SeqFeature' not in str(type(sub_feature)):
                        continue
                    else:
                        raw_features.append(sub_feature)
            else:
                continue
        prokka_rec.features = []
        added_cds = {}
        for feature in raw_features:
#            if feature.type == 'CDS':
#                check_locus = feature.qualifiers['locus_tag'][0]
#                if check_locus not in added_from_ratt and check_locus not in prokka_blast_list:
#                    print(check_locus)
            if feature.type != 'CDS':
                prokka_rec.features.append(feature)
            elif feature.type == 'CDS' and 'note' in feature.qualifiers.keys():
                if 'merged in this isolate' in feature.qualifiers['note'][0]:
                    if 'gene' in feature.qualifiers.keys():
                        feature.qualifiers.pop('gene')
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    added_cds[loc_key] = feature
                    prokka_rec.features.append(feature)
            if feature.type == 'CDS' and feature.qualifiers['locus_tag'][0][0] == 'L' and \
                    'gene' in feature.qualifiers.keys():
                if '_' in feature.qualifiers['gene'][0]:
                    counter += 1
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    if loc_key in added_cds.keys():
                        continue
                    feature_noref = prokka_noref_dict[loc_key]
#                    print(feature_noref)
                    if 'gene' in feature_noref.qualifiers.keys() and '_' in feature_noref.qualifiers['gene'][0]:
                        new_locus_tag = 'L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['locus_tag'] = [new_locus_tag]
                        new_protein_id = 'C:L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['protein_id'] = [new_protein_id]
                    elif 'gene' in feature_noref.qualifiers.keys() and '_' not in feature_noref.qualifiers['gene'][0]:
                        feature_noref.qualifiers['locus_tag'] = feature_noref.qualifiers['gene']
                        feature_noref.qualifiers.pop('gene')
                    else:
                        new_locus_tag = 'L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['locus_tag'] = [new_locus_tag]
                        new_protein_id = 'C:L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['protein_id'] = [new_protein_id]
                    added_cds[loc_key] = feature_noref
                    prokka_rec.features.append(feature_noref)
                else:
                    #print(feature)
#                    feature.qualifiers['locus_tag'] = feature.qualifiers['gene']
#                    feature.qualifiers.pop('gene')
                    #print(feature)
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    added_cds[loc_key] = feature
                    prokka_rec.features.append(feature)
            elif feature.type == 'CDS':
                loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                if loc_key not in added_cds.keys():
                    added_cds[loc_key] = feature
                    prokka_rec.features.append(feature)
                elif added_cds[loc_key].qualifiers['locus_tag'] == feature.qualifiers['locus_tag']:
                    continue
                else:
                    prokka_rec.features.append(feature)
#        ###############################################################################################
#        ####################### Blasting unannotated CDSs from Prokka with H37Rv ######################
#        ###############################################################################################
        mtb_fasta = os.environ['GROUPHOME'] + '/data/depot/hypotheticome-clinical/mtb-novel-genes.fasta'
        # Inverting key-value pair in dictionary
        rv_to_gene_name = dict((rv, gene_name) for gene_name, rv in gene_name_to_rv.iteritems())
        for feature in prokka_rec.features:
            if feature.type == 'CDS' and (('gene' not in feature.qualifiers.keys() and
                                          ('L_' in feature.qualifiers['locus_tag'][0] or
                                          'L2_' in feature.qualifiers['locus_tag'][0])) or
                                          ('gene' in feature.qualifiers.keys() and
                                           '_' in feature.qualifiers['gene'][0] and
                                           'L2_' in feature.qualifiers['locus_tag'][0])):
                query_sequence = feature.qualifiers['translation'][0]
                blast_to_mtb = NcbiblastpCommandline(subject=mtb_fasta, outfmt='"7 qseqid qlen sseqid slen qlen length'
                                                                               ' pident qcovs"')
                stdout_2, stderr_2 = blast_to_mtb(stdin=query_sequence)
                mtb_hit, all_hits_mtb = identify_top_hits(stdout_2, mtb=True)
                for gene in all_hits_mtb.keys():
                    note = 'locus_tag:' + gene + ':' + str(all_hits_mtb[gene])
                    if 'note' in feature.qualifiers.keys():
                        feature.qualifiers['note'].append(note)
                    else:
                        feature.qualifiers['note'] = [note]

        ###############################################################################################
        ####################### Adding translated sequence to RATT annotations  #######################
        ###############################################################################################
        for feature in prokka_rec.features:
            if feature.type == 'CDS' and 'translation' not in feature.qualifiers.keys():
                feature_sequence = translate(feature.extract(record_sequence), to_stop=True)
                feature.qualifiers['translation'] = [feature_sequence]

        ###############################################################################################
        ########################## Verifying annotations for essential genes  #########################
        ###############################################################################################
        essential_genes = get_essential_genes()
        all_rv_genes_in_isolate = []
        essential_genes_length = get_lengths_of_genes(essential_genes, h37rv_protein_fasta)
        annomerge_cds = 0
        cds_lengths = []
        for feature in prokka_rec.features:
            if feature.type == 'CDS':
                annomerge_cds += 1
                locus_tag = feature.qualifiers['locus_tag'][0]
                if locus_tag[:2] == 'Rv':
                    all_rv_genes_in_isolate.append(locus_tag)
                    if 'gene' not in feature.qualifiers.keys() and locus_tag in rv_to_gene_name.keys():
                        feature.qualifiers['gene'] = [rv_to_gene_name[locus_tag]]
                new_locus_tag = feature.qualifiers['locus_tag'][0]
                if 'gene' not in feature.qualifiers.keys():
                    feature.qualifiers['gene'] = [new_locus_tag]
                if new_locus_tag in essential_genes:
                    rv_length = essential_genes_length[new_locus_tag]
                    feature_length = len(feature.qualifiers['translation'][0])
                    feature_coverage = float(feature_length)/float(rv_length)
                    if feature_coverage <= 0.95:
                        print('WARNING: Essential gene ' + new_locus_tag + ' is truncated in isolate ' + args.isolate
                              + '. Coverage: ' + str(feature_coverage*100))
                cds_lengths.append(feature_length)
        output_file.write('Contig number: ' + str(i + 1) + '\n')
        output_file.write('Number of CDSs in isolate: ' + str(annomerge_cds) + '\n')
        output_file.write('Number of CDSs transferred from RATT: ' + str(cds_from_ratt) + '\n')
        output_file.write('Number of CDSs transferred from Prokka: ' + str(annomerge_cds-cds_from_ratt) + '\n')

        if len(cds_lengths) != 0:
            output_file.write(str('Min length of CDSs: ' + str(min(cds_lengths)) + '\n'))
            output_file.write(str('Max length of CDSs: ' + str(max(cds_lengths)) + '\n'))
            output_file.write(str('Median length of CDSs: ' + str(median(cds_lengths)) + '\n'))
        annomerge_records_post_processed.append(prokka_rec)
    SeqIO.write(annomerge_records_post_processed, output_genbank, 'genbank')


if __name__ == "__main__":
    main()
