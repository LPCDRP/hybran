#!/usr/bin/env python2.7
# Author: Deepika Gunasekaran
# Title: Merge annotation from Prokka for the positions which ar not annotated by 
# RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files
# in case of multiple contigs/chromosome annotations and a Genbank file (.gbf) file from Prokka annotation. 
# The output is an EMBL file with annotation predominantly from RATT and the intergenic regions annotated
# by RATT are filled with Prokka. This script also generates a log file to indicate characteristics of the 
# transferred features from RATT  


import sys
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import collections
from numpy import median
import os


def rename_locus(gene, strand):
    # This function checks names of existing locus tags in H37Rv and names the newly merged gene ensuring that the new
    # name does not conflict with existing locus tags
    char_start = 65 # ASCII conversion of 'A'
    gene_name = gene[:6]
    #print(gene)
    #print(gene_name)
    #print(strand)
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


def generate_feature_dictionary(feature_list):
    # This function takes as input a list of features and returns a dictionary with the
    # key as a tuple of feature start and stop positions and the value as the feature.
    feature_dict = {}
    for feature in feature_list:
        feature_key = (int(feature.location.start), int(feature.location.end), 
                       int(feature.location.strand))
        feature_dict[feature_key] = feature
    sorted_feature_dict = collections.OrderedDict(sorted(feature_dict.items()))
    return sorted_feature_dict


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


def check_inclusion_criteria(annotation_mapping_dict, embl_file, ratt_annotation, prokka_annotation):
    ratt_feature_start = int(ratt_annotation.location.start)
    ratt_feature_end = int(ratt_annotation.location.end)
    prokka_feature_start = int(prokka_annotation.location.start)
    prokka_feature_end = int(prokka_annotation.location.end)
    included = False
    # Check if feature types are the same. If not add feature to EMBL record
    if ratt_annotation.type != prokka_annotation.type:
        embl_file.features.append(prokka_annotation)
        included = True
    # Check if gene names match and if they don't or if gene names are missing, keep both
    elif 'gene' in ratt_annotation.qualifiers.keys() and 'gene' in prokka_annotation.qualifiers.keys():
        if ratt_annotation.qualifiers['gene'] != prokka_annotation.qualifiers['gene']:
            embl_file.features.append(prokka_annotation)
            included = True
        # If gene names are the same and the lengths of the genes are comparable between RATT and Prokka annotation
        # (difference in length of less than/equal to 10 bps), the RATT annotation is prefered
        elif ratt_annotation.qualifiers['gene'] == prokka_annotation.qualifiers['gene'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            embl_file.features.append(prokka_annotation)
            included = True
        # If gene tag is missing and the product is not a hypothetical protein, check to see if the products are the
            # same between RATT and Prokka and if they are and if the lengths of the protein coding genes are comparable
            # between RATT and Prokka (difference in length is less than/equal to 10 bps), the RATT annotation is
            # prefered
    elif ('gene' not in ratt_annotation.qualifiers.keys() or 'gene' not in prokka_annotation.qualifiers.keys()) and \
            ('product' in ratt_annotation.qualifiers.keys() and 'product' in prokka_annotation.qualifiers.keys()):
        if ratt_annotation.qualifiers['product'] == 'hypothetical protein' or \
                        prokka_annotation.qualifiers['product'] == 'hypothetical protein':
            embl_file.features.append(prokka_annotation)
            included = True
        elif ratt_annotation.qualifiers['product'] == prokka_annotation.qualifiers['product'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            embl_file.features.append(prokka_annotation)
            included = True
    return embl_file, included


def check_for_dnaA(feature_list):
    # This function takes as input a list of features and checks if the genome was circularized
    feature_dictionary = generate_feature_dictionary(feature_list)
    for cds in feature_dictionary.keys():
        if feature_dictionary[cds].type != 'CDS':
            continue
        if feature_dictionary[cds].qualifiers['locus_tag'][0] != 'Rv0001':
            #print(feature_dictionary[cds])
            print('DnaA is not the first gene in this genome. Possible circularization error')
            sys.exit()
        elif feature_dictionary[cds].qualifiers['locus_tag'][0] == 'Rv0001' and \
                int(feature_dictionary[cds].location.start) == 0:
            #print(int(feature_dictionary[cds].location.start))
            break
        else:
            print('DnaA is the first gene in this genome. But the position is off')
            sys.exit()
    return


def isolate_valid_ratt_annotations(feature_list):
    # This function takes as input a list of features and checks if the length of the CDSs are divisible by 3. If not,
    # it outputs the features to stdout and removes them from the valid_ratt_annotations
    valid_features = []
    for feature in feature_list:
        if feature.type == 'CDS' and (len(feature.location) % 3) != 0:
            print(feature)
        else:
            valid_features.append(feature)
    return valid_features



# Required inputs:
# 1. Isolate ID or Explicit RATT annotation EMBL file(s) and Prokka annotation Genbank file
# Optional inputs:
# 1. Output file/path for Genbank output
# 2. Output file/path for log file


def main():
    input_ratt_files = []
    # Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file>
    # from the commandline. The log file output stats about the features that are added to the RATT 
    # annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/
    # 'isolate_id'.gbf and 'isolate_id'/annomerge/'isolate_id'.log
    parser = argparse.ArgumentParser(description='Merging annotation from RATT and Prokka',
                                     epilog='Isolate ID must be specified OR Explicit file paths for RATT and Prokka '
                                            'annotation must be specified')
    parser.add_argument('-i', '--isolate', help='Isolate ID')
    parser.add_argument('-o', '--output', help='Output file in Genbank format', default='annomerge.gbf')
    parser.add_argument('-l', '--log_file', help='Log file with information on features added from prokka',
                        default='annomerge_annotation.log')
    parser.add_argument('-m', '--merged_genes', help='Merged genes file in genbank format',
                        default='merged_genes.gbf')
    parser.add_argument('-r', '--ratt_annotation', nargs='*', help='RATT annotation file(s) (i.e. .final.embl) file(s)')
    parser.add_argument('-p', '--prokka_annotation', help='Prokka annotation Genbank file')
    args = parser.parse_args()
    if args.isolate is None and (args.ratt_annotation is None or args.prokka_annotation is None):
        parser.print_help()
        sys.exit('INVALID INPUT')
    if args.isolate is not None:
        file_path = os.environ['GROUPHOME'] + '/data/depot/annotation/' + args.isolate + '/'
        ratt_file_path = file_path + 'ratt'
        try:
            ratt_embl_files = [embl_file for embl_file in os.listdir(ratt_file_path)
                               if embl_file.endswith('.final.embl')]
            for embl_file in ratt_embl_files:
                embl_file_path = ratt_file_path + '/' + embl_file
                input_ratt_files.append(embl_file_path)
        except OSError:
            sys.exit('Expecting RATT annotation files but found none')
        try:
            input_prokka_genbank = file_path + 'prokka/' + args.isolate + '.gbf'
        except OSError:
            sys.exit('Expecting Prokka annotation file but found none')
    if args.isolate is None:
        input_ratt_files = args.ratt_annotation
        input_prokka_file = args.prokka_annotation
        if not os.path.isfile(input_prokka_file):
            sys.exit('Invalid file path for Prokka annotation')
        for ratt_input in input_ratt_files:
            if not os.path.isfile(ratt_input):
                sys.exit('Invalid file path for RATT annotation')
    output_merged_genes = args.merged_genes
    output_genbank = args.output
    output_log =  args.log_file
    output_file = open(output_log, 'w+')
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    annomerge_records = []
    global genes_in_rv
    genes_in_rv = []
    genes_in_rv_file = os.environ['GROUPHOME'] + '/resources/H37rv_proteincoding_genes.txt'
    genes_in_rv_raw = open(genes_in_rv_file, 'r').readlines()
    for line in genes_in_rv_raw:
        locus_tag = line.split()[0]
        genes_in_rv.append(locus_tag)
    for i in range(0, len(input_ratt_files)):
        ratt_contig_record = SeqIO.read(input_ratt_files[i], 'embl')
        prokka_contig_record = prokka_records[i]
        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        if i == 0:
            check_for_dnaA(ratt_contig_features)
        ratt_contig_features = isolate_valid_ratt_annotations(ratt_contig_features)
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
        if len(merged_features) > 0 and i ==0:
            SeqIO.write(merged_features_record, output_merged_genes, 'genbank')
        if len(ratt_contig_features) == 0:
            print("NO RATT ANNOTATION FOR CONTIG " + str(i+1))
            feature_additions = {}
            feature_lengths = {}
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
            #print("Annomerge record length")
            #print(len(annomerge_records))
            for prokka_feature in prokka_contig_record.features:
                if prokka_feature.type not in feature_additions.keys():
                    feature_additions[prokka_feature.type] = 1
                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                else:
                    feature_additions[prokka_feature.type] += 1
                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            output_file.write('Contig Number: ' + str(i+1) + '\n')
            for f in feature_lengths.keys():
                output_file.write(str('Feature: ' + f + '\n'))
                output_file.write(str('Number of ' + f + ': ' + str(feature_additions[f]) + '\n'))
                output_file.write(str('Min length: ' + str(min(feature_lengths[f])) + '\n'))
                output_file.write(str('Max length: ' + str(max(feature_lengths[f])) + '\n'))
                output_file.write(str('Median length: ' + str(median(feature_lengths[f])) + '\n'))
            #annomerge_contig_record = prokka_contig_record
            continue
        elif len(prokka_contig_features) == 0:
            print("NO PROKKA ANNOTATION FOR CONTIG " + str(i+1))
            #annomerge_contig_record = ratt_contig_record
            prokka_contig_features = ratt_contig_features
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            output_file.write('Contig Number: ' + str(i+1) + '\n')
            output_file.write('No Annotation to add from Prokka')
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the
            # features from RATT
            annomerge_contig_record = prokka_contig_record[:]
            annomerge_contig_record.features = ratt_contig_features
            annomerge_contig_record.annotations['comment'] = 'Merged reference based annotation from RATT and ab ' \
                                                             'initio annotation from Prokka'
            ###########################################################
            ####### Creating a dictionary with feature location #######
            ############### as key and index as value #################
            ###########################################################
            ratt_annotation_mapping = {} # Used for resolving annotations of overlapping features between RATT and
            # Prokka
            for index, feature in enumerate(ratt_contig_record.features):
                start = int(feature.location.start)
                end = int(feature.location.end)
                ratt_annotation_mapping[(start, end)] = index
            prokka_annotation_mapping = {}
            ratt_contig_record_mod = ratt_contig_record[:]
            ratt_contig_record_mod.features = ratt_contig_features
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
                            prokka_features_not_in_ratt.pop((prokka_feature_start,prokka_feature_end,prokka_strand),
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
                            annomerge_contig_record.features.append(prokka_feature)
                        # If the Prokka feature overlaps with two RATT features
                        elif prokka_feature_start < ratt_unannotated_region_start and \
                                        prokka_feature_end > ratt_unannotated_region_end:
                            if prokka_feature.type == 'source': # This is to exclude the source feature as it is
                                # accounted for by RATT
                                break
                            annomerge_contig_record.features.append(prokka_feature)
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
                                ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start-1,
                                                                               prokka_strand)]
                            elif (prokka_feature_start in ratt_unannotated_region_range) and \
                                            prokka_feature_end > ratt_unannotated_region_end:
                                ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                                prokka_strand)]
                            annomerge_contig_record, included = check_inclusion_criteria(ratt_annotation_mapping,
                                                                                         annomerge_contig_record,
                                                                                         ratt_overlapping_feature,
                                                                                         prokka_feature)
                            if included: # To check if Prokka feature passed the inclusion criteria and was integrated
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
                for feature in merged_features:
                    annomerge_contig_record.features.append(feature)
            annomerge_records.append(annomerge_contig_record)
            output_file.write('Contig Number: ' + str(i+1) + '\n')
            for f in feature_lengths.keys():
                output_file.write(str('Feature: ' + f + '\n'))
                output_file.write(str('Number of ' + f + ': ' + str(feature_additions[f]) + '\n'))
                output_file.write(str('Min length: ' + str(min(feature_lengths[f])) + '\n'))
                output_file.write(str('Max length: ' + str(max(feature_lengths[f])) + '\n'))
                output_file.write(str('Median length: ' + str(median(feature_lengths[f])) + '\n'))
    output_file.write(str('Number of merged genes: ' + str(len(merged_features)) + '\n'))
    SeqIO.write(annomerge_records, output_genbank, 'genbank')

if __name__ == "__main__":
   main()
