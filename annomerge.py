#!/usr/bin/env python2.7
# Author: Deepika Gunasekaran
# Title: Merge annotation from Prokka for the positions which ar not annotated by 
# RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation
# and a Genbank file (.gbf) file from Prokka annotation. The output is an EMBL file
# with annotation predominantly from RATT and the intergenic regions annotated by 
# RATT are filled with Prokka. This script also generates a log file to indicate
# characteristics of the transferred features from RATT  


import sys, getopt
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from numpy import median
import collections

#################################################################################
# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
# This function was modified by Deepika Gunasekaran
def get_interregions(genbank_path,intergene_length=1):
    seq_record = SeqIO.parse(open(genbank_path), "embl").next()
    cds_list_plus = []
    cds_list_minus = []
    intergenic_records = []
    intergenic_positions = []
    # Loop over the genome file, get the CDS features on each of the strands
    for feature in seq_record.features:
#        if feature.type == 'CDS':
        mystart = feature.location.start.position
        myend = feature.location.end.position
        if feature.strand == -1:
            cds_list_minus.append((mystart,myend,-1))
        elif feature.strand == 1:
            cds_list_plus.append((mystart,myend,1))
        else:
            sys.stderr.write("No strand indicated %d-%d. Assuming +\n" %
                                  (mystart, myend))
            cds_list_plus.append((mystart,myend,1))
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
    return intergenic_records, intergenic_positions
#####################################################################################


def generate_feature_dictionary(feature_list):
# This function takes as input a list of features and returns a dictionary with the key as a tuple of
# feature start and stop positions and the value as the feature.
    feature_dict = {}
    for feature in feature_list:
        feature_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
        feature_dict[feature_key] = feature
    sorted_feature_dict = collections.OrderedDict(sorted(feature_dict.items()))
    return sorted_feature_dict


def remove_duplicate_annotations(ratt_features, prokka_features_dictionary):
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




# Required input files:
# 1. RATT annotation in a valid embl file.
# 2. Prokka annotation in a gff file.

def main(argv):
    input_ratt_embl = ''
    input_prokka_genbank = ''
    output_embl = ''
    #output_gff = ''
    output_log = ''
    # Annomerge takes as options -r <input_ratt_embl_file> -p <input_prokka_genbank_file> -o <output_file>
    # from the commandline. Additionally, it also takes as option -l <ouput_log_file>
    # in case a lof file is required for troubleshooting. The log file output stats about the features that
    # are added to the RATT annotation
    try:
        opts, args = getopt.getopt(argv,"hr:p:e:l:",["ratt=", "prokka=", "output_embl=", "log_file="])
    except getopt.GetoptError:
        print 'annomerge.py -r <input_ratt_embl_file> -p <input_prokka_genbank_file> -e <output_embl_file> -l <output_log_file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'annomerge.py -r <input_ratt_embl_file> -p <input_prokka_genbank_file> -e <output_embl_file> -l <output_log_file>'
            sys.exit()
        elif opt in ("-r", "--ratt"):
            input_ratt_embl = arg
        elif opt in ("-p", "--prokka"):
            input_prokka_genbank = arg
        elif opt in ("-e", "--output_embl"):
            output_embl = arg
        elif opt in ("-l", "--log_file"):
            output_log = arg

    ratt_record = SeqIO.read(input_ratt_embl, 'embl')
    #print('RATT record: ' + str(len(ratt_record.features)))
    prokka_record = SeqIO.read(input_prokka_genbank, 'genbank')
    #print('Prokka record: ' + str(len(prokka_record.features)))
    embl_record = SeqIO.read(input_ratt_embl, 'embl')
    #print(len(embl_record.features))
    #print('RATT annotation')
    ratt_cds_pos = {}
    ratt_cds_neg = {}
    prokka_cds_pos = {}
    prokka_cds_neg = {}
    ratt_gaps = []
    prokka = {}
    ratt_features = ratt_record.features
    prokka_features = prokka_record.features
    prokka_features_dict = generate_feature_dictionary(prokka_features)
    #print('Original Prokka feature length: ' + str(len(prokka_features_dict.keys())))
    prokka_features_not_in_ratt, ratt_overlapping_genes = remove_duplicate_annotations(ratt_features, prokka_features_dict)
    #print('Reduced Prokka feature length: ' + str(len(prokka_features_not_in_ratt.keys())))
    #print('Potential RATT gene overlaps: ' + str(len(ratt_overlapping_genes.keys())))
    intergenic_ratt, intergenic_positions = get_interregions(input_ratt_embl, intergene_length=1)
    feature_additions = {}
    feature_lengths = {}
    overlap_features = {}
    overlap_feature_lengths = {}
    #overlap_feature_overlap_percent = []
    corner_cases = 0
    for i in intergenic_positions:
        # Variable definitions
        ratt_start = i[0]
        ratt_end = i[1]
        ratt_strand = i[2]
        for feature_position in prokka_features_dict.keys():
            prokka_strand = feature_position[2]
            prokka_start = feature_position[0]
            prokka_end = feature_position[1]
            ratt_pos_range = range(ratt_start, ratt_end + 1)
            if (prokka_strand == -1 and ratt_strand == '-') or (prokka_strand == 1 and ratt_strand == '+'):
                # If Prokka feature is location before the start of the intergenic region, pop the key from the dictionary
                if (prokka_start < ratt_start) and (prokka_end < ratt_start):
                    prokka_features_dict.pop((prokka_start,prokka_end,prokka_strand), None)
                    continue
                # Else if the prokka feature location is after the end of the intergenic region, break out of the inner loop
                elif (prokka_start > ratt_end) and (prokka_end > ratt_end):
                    break
                # If the PROKKA feature is contained in the RATT feature
                elif prokka_start in ratt_pos_range and prokka_end in ratt_pos_range:
                    feature = prokka_features_dict[feature_position]
                    if feature.type not in feature_additions.keys():
                        feature_additions[feature.type] = 1
                        feature_lengths[feature.type] = [len(feature.location)]
                    else:
                        feature_additions[feature.type] += 1
                        feature_lengths[feature.type].append(len(feature.location))
                    embl_record.features.append(feature)

               #elif int(prokka_start) < ratt_start and int(prokka_end) in ratt_pos_range:
                    #if feature.type not in overlap_features.keys():
                        #overlap_features[feature.type] = 1
                        #overlap_feature_lengths[feature.type] = [len(feature.location)]
                    #else:
                        #overlap_features[feature.type] += 1
                        #overlap_feature_lengths[feature.type].append(len(feature.location))
                    #overlap_percent = (float(len(feature.location)) - float(ratt_start)) / float(len(feature.location))
                    #overlap_feature_overlap_percent.append(overlap_percent)
                    #p_overlap = (ratt_start - int(prokka_start) / (len(feature.location)))
                #elif int(prokka_start) in ratt_pos_range and int(prokka_end) > ratt_end:
                    #if feature.type not in overlap_features.keys():
                        #overlap_features[feature.type] = 1
                        #overlap_feature_lengths[feature.type] = [len(feature.location)]
                    #else:
                        #overlap_features[feature.type] += 1
                        #overlap_feature_lengths[feature.type].append(len(feature.location))
                    #overlap_percent = (float(len(feature.location))-ratt_end)/float(len(feature.location))
                    #overlap_feature_overlap_percent.append(overlap_percent)
                    #p_overlap = (int(prokka_start)-ratt_end)/(len(feature.location))
                #elif int(prokka_start) < ratt_start and int(prokka_end) > ratt_end:
                    #corner_cases += 1
            # else:
            #    continue
    #print(feature_additions)
    if output_log:
        output_file = open(output_log, 'w')
        for f in feature_lengths.keys():
            output_file.write(str('Feature: ' + f + '\n'))
            output_file.write(str('Number of ' + f + ': ' + str(feature_additions[f]) + '\n'))
            output_file.write(str('Min length: ' + str(min(feature_lengths[f])) + '\n'))
            output_file.write(str('Max length: ' + str(max(feature_lengths[f])) + '\n'))
            output_file.write(str('Median length: ' + str(median(feature_lengths[f])) + '\n'))
    #print(len(embl_record.features))
    if output_embl:
        print('Writing Output EMBL file')
        SeqIO.write(embl_record, output_embl, 'embl')
    #for f in overlap_feature_lengths.keys():
        #print(str('Feature: ' + f + '\n'))
        #print(str('Number of ' + f + ': ' + str(overlap_features[f]) + '\n'))
        #print(str('Min length: ' + str(min(overlap_feature_lengths[f])) + '\n'))
        #print(str('Max length: ' + str(max(overlap_feature_lengths[f])) + '\n'))
        #print(str('Median length: ' + str(median(overlap_feature_lengths[f])) + '\n'))
    #print(str('Corner cases: ' + str(corner_cases) + '\n'))
    #for i in overlap_feature_overlap_percent:
        #print(i)


if __name__ == "__main__":
   main(sys.argv[1:])
  
