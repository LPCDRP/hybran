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

    ratt_record = SeqIO.parse(input_ratt_embl, 'embl')
    prokka_record = SeqIO.parse(input_prokka_genbank, 'genbank')
    embl_record = SeqIO.read(input_ratt_embl, 'embl')
    #print(len(embl_record.features))
    #print('RATT annotation')
    ratt_cds_pos = {}
    ratt_cds_neg = {}
    prokka_cds_pos = {}
    prokka_cds_neg = {}
    ratt_gaps = []
    prokka = {}
    for rec in ratt_record:
        ratt_features = rec.features
    for rec in prokka_record:
        prokka_features = rec.features
    intergenic_ratt,intergenic_positions = get_interregions(input_ratt_embl,intergene_length=1)
    feature_additions = {}
    feature_lengths = {}
    overlap_features = {}
    overlap_feature_lengths = {}
    #overlap_feature_overlap_percent = []
    corner_cases = 0

    for feature in prokka_features:
        start = feature.location.start
        stop = feature.location.end
        key = (start, stop)
        prokka[key] = feature
    for i in intergenic_positions:
        #print(i[2])
        ratt_strand = i[2]
        feature = prokka[(i[0], i[1])]
        # for feature in prokka_features:
        prokka_strand = feature.strand
        if prokka_strand == -1 and ratt_strand == '-':
            if int(feature.location.start) in range(i[0],i[1]+1) and int(feature.location.end) in \
                    range(i[0],i[1]+1):
                #print(feature)
                if feature.type not in feature_additions.keys():
                    feature_additions[feature.type] = 1
                    feature_lengths[feature.type] = [len(feature.location)]
                else:
                    feature_additions[feature.type] += 1
                    feature_lengths[feature.type].append(len(feature.location))
                #print(len(feature.location))
                embl_record.features.append(feature)
            elif int(feature.location.start) < i[0] and int(feature.location.end) in range(i[0],i[1]+1):
                if feature.type not in overlap_features.keys():
                    overlap_features[feature.type] = 1
                    overlap_feature_lengths[feature.type] = [len(feature.location)]
                else:
                    overlap_features[feature.type] += 1
                    overlap_feature_lengths[feature.type].append(len(feature.location))
                #overlap_feature_overlap_percent.append((float(len(feature.location))-float(i[0]))/float(len(feature.location)))
                p_overlap = (i[0]-int(feature.location.start)/(len(feature.location)))
                #print(p_overlap)
                #print(feature)
            elif int(feature.location.start) in range(i[0],i[1]+1) and int(feature.location.end) > i[1]:
                if feature.type not in overlap_features.keys():
                    overlap_features[feature.type] = 1
                    overlap_feature_lengths[feature.type] = [len(feature.location)]
                else:
                    overlap_features[feature.type] += 1
                    overlap_feature_lengths[feature.type].append(len(feature.location))
                #print(len(feature.location))
                #overlap_feature_overlap_percent.append((float(len(feature.location))-i[1])/float(len(feature.location)))
                p_overlap = (int(feature.location.start)-i[1])/(len(feature.location))
                #print(p_overlap)
                #print(feature)
            elif int(feature.location.start) < i[0] and int(feature.location.end) > i[1]:
                corner_cases += 1
        elif prokka_strand == 1 and ratt_strand == '+':
            if int(feature.location.start) in range(i[0],i[1]+1) and int(feature.location.end) in range(i[0],i[1]+1):
                #print(feature)
                if feature.type not in feature_additions.keys():
                    feature_additions[feature.type] = 1
                    feature_lengths[feature.type] = [len(feature.location)]
                else:
                    feature_additions[feature.type] += 1
                    feature_lengths[feature.type].append(len(feature.location))
                embl_record.features.append(feature)
            elif int(feature.location.start) < i[0] and int(feature.location.end) in range(i[0],i[1]+1):
                if feature.type not in overlap_features.keys():
                    overlap_features[feature.type] = 1
                    overlap_feature_lengths[feature.type] = [len(feature.location)]
                else:
                    overlap_features[feature.type] += 1
                    overlap_feature_lengths[feature.type].append(len(feature.location))
                #overlap_feature_overlap_percent.append((float(len(feature.location))-i[0])/float(len(feature.location)))
                p_overlap = (i[0]-int(feature.location.start)/(len(feature.location)))
                #print(p_overlap)
                #print(feature)
            elif int(feature.location.start) in range(i[0],i[1]+1) and int(feature.location.end) > i[1]:
                if feature.type not in overlap_features.keys():
                    overlap_features[feature.type] = 1
                    overlap_feature_lengths[feature.type] = [len(feature.location)]
                else:
                    overlap_features[feature.type] += 1
                    overlap_feature_lengths[feature.type].append(len(feature.location))
                #print(len(feature.location))
                #overlap_feature_overlap_percent.append((float(len(feature.location))-i[1])/float(len(feature.location)))
                p_overlap = (int(feature.location.start)-i[1])/(len(feature.location))
                #print(p_overlap)
                #print(feature)
            elif int(feature.location.start) < i[0] and int(feature.location.end) > i[1]:
                corner_cases += 1
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
  
