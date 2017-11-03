# Author: Deepika Gunasekaran
# Title: Convert EMBL file with nucleotide sequence to fasta file with protein sequence
# Description: This program takes as input, an EMBL file of H37Rv with mannotation and converts it into a# fasta file with protein sequences

import sys
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import os


def main():
    input_ratt_files = []
    # Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file>
    # from the commandline. The log file output stats about the features that are added to the RATT 
    # annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/
    # 'isolate_id'.gbf and 'isolate_id'/annomerge/'isolate_id'.log
    parser = argparse.ArgumentParser(description='Converting EMBL file with nucleotide sequence to FASTA file with amino acid sequence')
    parser.add_argument('-e', '--input_embl', help='Input EMBL file', required=True)
    parser.add_argument('-o', '--output_fasta', help='Output FASTA file', required=True)
    args = parser.parse_args()
    input_embl = args.input_embl
    output =  args.output_fasta
    output_file = open(output, 'wb')
    embl_record = SeqIO.read(input_embl, 'embl')
    id_entry = embl_record.id
    for feature in embl_record.features:
        if feature.type == 'CDS':
            #fasta_header = '>gi|' + id_entry + '|gb|' + id_entry + '|' + feature.qualifiers['locus_tag'][0]
            if 'product' in feature.qualifiers.keys() and 'translation' in feature.qualifiers.keys():
                product_list = ' '.join(feature.qualifiers['product'])
                if 'gene' in feature.qualifiers.keys():
                    fasta_header = '>1 ' + '~~~' + feature.qualifiers['gene'][0] + '~~~' + product_list
                    gene_list = ' '.join(feature.qualifiers['gene'])
                    fasta_header = fasta_header + ' ' + gene_list + ' ' + product_list + '\n'
                else:
                    fasta_header = '>1 ' + '~~~' + feature.qualifiers['locus_tag'][0] + '~~~' + product_list
                    fasta_header = fasta_header + ' ' + product_list + '\n'
                output_file.write(fasta_header)
                output_file.write(feature.qualifiers['translation'][0])
                output_file.write('\n')
            #else:
                #print(feature)


if __name__ == "__main__":
   main()
