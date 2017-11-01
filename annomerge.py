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
        if feature.type == 'mRNA':
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
        if ratt_feature.type == 'gene':
            continue
        ratt_start = ratt_feature.location.start
        ratt_end = ratt_feature.location.end
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
                if ('gene' in ratt_feature.qualifiers.keys() and 'gene' in prokka_feature.qualifiers.keys()) and (ratt_feature.qualifiers['gene'] == prokka_feature.qualifiers['gene']):
                    prokka_features_not_in_ratt.pop((prokka_start,prokka_end,prokka_strand), None)
                    prokka_duplicate_removed = True
                elif len(ratt_feature.location) == len(prokka_feature.location):
                    prokka_features_not_in_ratt.pop((prokka_start,prokka_end,prokka_strand), None)
                    prokka_duplicate_removed = True
        if not prokka_duplicate_removed and ratt_feature.type != 'mRNA':
            ratt_overlapping_genes[(ratt_start,ratt_end,ratt_strand)] = ratt_feature
    return prokka_features_not_in_ratt, ratt_overlapping_genes


def check_inclusion_criteria(annotation_mapping_dict, embl_file, ratt_annotation, prokka_annotation):
    ratt_feature_start = ratt_annotation.location.start
    ratt_feature_end = ratt_annotation.location.end
    prokka_feature_start = prokka_annotation.location.start
    prokka_feature_end = prokka_annotation.location.end
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
        # If gene names are the same and the lengths of the genes are comparable between RATT and Prokka annotation (difference in length of less than/equal to 10 bps), the RATT annotation is prefered
        elif ratt_annotation.qualifiers['gene'] == prokka_annotation.qualifiers['gene'] and abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            embl_file.features.append(prokka_annotation)
            included = True
        # If gene tag is missing and the product is not a hypothetical protein, check to see if the products are the same between RATT and Prokka and if they are and if the lengths of the protein coding genes are comparable between RATT and Prokka (difference in length is less than/equal to 10 bps), the RATT annotation is prefered
    elif ('gene' not in ratt_annotation.qualifiers.keys() or 'gene' not in prokka_annotation.qualifiers.keys()) and ('product' in ratt_annotation.qualifiers.keys() and 'product' in prokka_annotation.qualifiers.keys()):
        if ratt_annotation.qualifiers['product'] == 'hypothetical protein' or prokka_annotation.qualifiers['product'] == 'hypothetical protein':
            embl_file.features.append(prokka_annotation)
            included = True
        elif ratt_annotation.qualifiers['product'] == prokka_annotation.qualifiers['product'] and abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 10:
            embl_file.features.append(prokka_annotation)
            included = True
    return embl_file, included
    


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
    parser = argparse.ArgumentParser(description='Merging annotation from RATT and Prokka', epilog='Isolate ID must be specified OR Explicit file paths for RATT and Prokka annotation must be specified')
    parser.add_argument('-i', '--isolate', help='Isolate ID')
    parser.add_argument('-o', '--output', help='Output file in Genbank format', default='annomerge.gbf')
    parser.add_argument('-l', '--log_file', help='Log file with information on features added from prokka', default='annomerge.log')
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
            ratt_embl_files = [embl_file for embl_file in os.listdir(ratt_file_path) if embl_file.endswith('.final.embl')]
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
    output_genbank = args.output
    output_log =  args.log_file
    output_file = open(output_log, 'w+')
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    annomerge_records = []
    for i in range(0,len(input_ratt_files)):
        ratt_contig_record = SeqIO.read(input_ratt_files[i], 'embl')
        prokka_contig_record = prokka_records[i]
        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        if len(ratt_contig_features) == 0:
            print("NO RATT ANNOTATION")
            feature_additions = {}
            feature_lengths = {}
            annomerge_records.append(prokka_contig_record)
            print("Annomerge record length")
            print(len(annomerge_records))
            for prokka_feature in prokka_contig_features:
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
            annomerge_contig_record = prokka_contig_record
            continue
        elif len(prokka_contig_features) == 0:
            print("NO PROKKA ANNOTATION")
            annomerge_contig_record = ratt_contig_record
            prokka_contig_features = ratt_contig_features
            annomerge_records.append(prokka_contig_record)
            print(len(prokka_contig_feature))
            output_file.write('Contig Number: ' + str(i+1) + '\n')
            output_file.write('No Annotation to add from Prokka')            
            for ratt_contig_feature in ratt_contig_features:
                prokka_contig_features.append(ratt_contig_feature)
            print(len(ratt_contig_features))
            print(len(annomerge_records))
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the features from RATT
            annomerge_contig_record = prokka_contig_record
            annomerge_contig_record.features = ratt_contig_features
            annomerge_contig_record.annotations['comment'] = 'Merged reference based annotation from RATT and ab initio annotation from Prokka'

            ###########################################################
            ####### Creating a dictionary with feature location #######
            ############### as key and index as value #################
            ###########################################################
            ratt_annotation_mapping = {} # Used for resolving annotations of overlapping features between RATT and Prokka
            for index, feature in enumerate(ratt_contig_record.features):
                start = feature.location.start
                end = feature.location.end
                ratt_annotation_mapping[(start,end)] = index
            prokka_annotation_mapping = {}
    
            prokka_features_dict = generate_feature_dictionary(prokka_contig_features)
            prokka_features_not_in_ratt, ratt_overlapping_genes = remove_duplicate_annotations(ratt_contig_features, prokka_features_dict)
            intergenic_ratt, intergenic_positions, ratt_pre_intergene, ratt_post_intergene = get_interregions(ratt_contig_record, intergene_length=1)
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
                    ratt_unannotated_region_range = range(ratt_unannotated_region_start, ratt_unannotated_region_end + 1)
                    if (prokka_strand == -1 and ratt_strand == '-') or (prokka_strand == 1 and ratt_strand == '+'):
                        # If Prokka feature is location before the start of the intergenic region, pop the key from the dictionary and continue loop
                        if (prokka_feature_start < ratt_unannotated_region_start) and (prokka_feature_end < ratt_unannotated_region_start):
                            prokka_features_not_in_ratt.pop((prokka_feature_start,prokka_feature_end,prokka_strand), None)
                            continue
                        # Else if the prokka feature location is after the end of the intergenic region, break out of the inner loop
                        elif (prokka_feature_start > ratt_unannotated_region_end) and (prokka_feature_end > ratt_unannotated_region_end):
                            break
                        # If the PROKKA feature is contained in the RATT feature
                        elif prokka_feature_start in ratt_unannotated_region_range and prokka_feature_end in ratt_unannotated_region_range:
                            # The if-else condition below is to keep track of the features added from Prokka for the log file
                            if prokka_feature.type not in feature_additions.keys():
                                feature_additions[prokka_feature.type] = 1
                                feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                            else:
                                feature_additions[prokka_feature.type] += 1
                                feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
                            annomerge_contig_record.features.append(prokka_feature)
                        # If the Prokka feature overlaps with two RATT features
                        elif prokka_feature_start < ratt_unannotated_region_start and prokka_feature_end > ratt_unannotated_region_end:
                            if prokka_feature.type == 'source': # This is to exclude the source feature as it is accounted for by RATT
                                break
                            annomerge_contig_record.features.append(prokka_feature)
                            # The if-else condition below is to keep track of the features added from Prokka for the log file
                            if prokka_feature.type not in feature_additions.keys():
                                feature_additions[prokka_feature.type] = 1
                                feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                            else:
                                feature_additions[prokka_feature.type] += 1
                                feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
                        # If the Prokka feature overlaps with one RATT feature
                        else:
                            if (prokka_feature_start < ratt_unannotated_region_start) and (prokka_feature_end in ratt_unannotated_region_range): 
                                ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start-1,prokka_strand)]
                            elif (prokka_feature_start in ratt_unannotated_region_range) and prokka_feature_end > ratt_unannotated_region_end:
                                ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,prokka_strand)]
                            annomerge_contig_record, included = check_inclusion_criteria(ratt_annotation_mapping, annomerge_contig_record, ratt_overlapping_feature, prokka_feature)
                            if included: # To check if Prokka feature passed the inclusion criteria and was integrated into the EMBL file
                            # The if-else condition below is to keep track of the features added from Prokka for the log file
                                if prokka_feature.type not in feature_additions.keys():
                                    feature_additions[prokka_feature.type] = 1
                                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                                else:
                                    feature_additions[prokka_feature.type] += 1
                                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            annomerge_records.append(annomerge_contig_record)

            output_file.write('Contig Number: ' + str(i+1) + '\n')
            for f in feature_lengths.keys():
                output_file.write(str('Feature: ' + f + '\n'))
                output_file.write(str('Number of ' + f + ': ' + str(feature_additions[f]) + '\n'))
                output_file.write(str('Min length: ' + str(min(feature_lengths[f])) + '\n'))
                output_file.write(str('Max length: ' + str(max(feature_lengths[f])) + '\n'))
                output_file.write(str('Median length: ' + str(median(feature_lengths[f])) + '\n'))

    SeqIO.write(annomerge_records, output_genbank, 'genbank')

if __name__ == "__main__":
   main()
  
