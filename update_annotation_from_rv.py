#!/usr/bin/env python2.7
# Author: Deepika Gunasekaran
# Title: Update annotation for novel genes from H37Rv
# Description: This program takes as input, the genbank file from annomerge outputs, compiles a list of CDS candidates
# that do not have an annotated gene name and performs a blastp with H37Rv protein fasta file. If the amino acid
# identity is greater than a certain threshold (default: 85%), this script updates the annotation of the genbank file
# with the corresponding H37Rv hit.

### TO PURGE IN TEMP DIRECTORY:
# 1. unannotated_genes.faa
# 2. blast_output_file

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import sys
import argparse
from argparse import RawTextHelpFormatter


def genesLocusTagParser(file_path=None):
    """
    :param file_path: Filepath for genes.tsv where 1st column is Rv locus tag and second column is gene name if present.
    Default file used is $GROUPHOME/resources/mtb-reconstruction/genes.tsv
    :return: Dictionary where keys are Rv numbers and values are gene names. Only Rvs with corresponding gene names are
    present in this dictionary
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
            genes_dict[line_elements[0]] = line_elements[1]
    return genes_dict


def getUnannotatedGenes(genbank_filename):
    """
    :param genbank_filename: Annotated Genbank file from Annomerge
    :return: Fasta file path for unannotated genes fasta file
    """
    fasta_filename = output_dir_name + '/unannotated_genes.faa'
    output_handle = open(fasta_filename, 'w')
    genbank_filepath = args.annotation_directory + '/' + genbank_filename
    for record in SeqIO.parse(genbank_filepath, 'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                if feature.qualifiers['locus_tag'][0][0] == 'L' and 'gene' not in feature.qualifiers.keys():
                    output_handle.write(">%s from %s\n%s\n" % (feature.qualifiers['locus_tag'][0], record.name,
                                                               feature.qualifiers['translation'][0]))
    output_handle.close()
    return fasta_filename


def filterTopHits(all_hits_dict):
    """
    :param all_hits_dict: This dictionary consists of all Blast hits with corresponding unannotated locus_tag that have
    passed the defined thresholds for amino acid identity and coverage
    :return: returns only the elements of top hits with the locus tag. The key is the 'L(2)_' locus tag and the value is
     corresponding top hit in H37Rv. If multiple top hits are present with same identity and coverage values, the value
     is all the top Rv hits separated by ':'
    """
    top_hits = {}
    for locus_tag in all_hits_dict.keys():
        if all_hits_dict[locus_tag] == {}:
            continue
        elif len(all_hits_dict[locus_tag].keys()) == 1:
            top_hits[locus_tag] = all_hits_dict[locus_tag].keys()[0]
        else:
            # Top hits with higher coverage considered over higher amino acid % identity
            max_cov_dict = {}
            max_cov = 0
            for match in all_hits_dict[locus_tag].keys():
                if all_hits_dict[locus_tag][match][1] > max_cov:
                    max_cov = all_hits_dict[locus_tag][match][1]
                    max_cov_dict = {match: all_hits_dict[locus_tag][match]}
                elif all_hits_dict[locus_tag][match][1] == max_cov:
                    max_cov_dict[match] = all_hits_dict[locus_tag][match]
                else:
                    continue
            if len(max_cov_dict.keys()) == 1:
                top_hits[locus_tag] = max_cov_dict.keys()[0]
            else:
                max_iden = 0
                max_iden_dict = {}
                for gene in max_cov_dict.keys():
                    if max_cov_dict[gene][0] > max_iden:
                        max_iden = max_cov_dict[gene][0]
                        max_iden_dict = {gene: max_cov_dict[gene]}
                    elif max_cov_dict[gene][0] == max_iden:
                        max_iden_dict[gene] = max_cov_dict[gene]
                    else:
                        continue
            if len(max_iden_dict.keys()) == 1:
                top_hits[locus_tag] = max_iden_dict.keys()[0]
            else:
                print('Multiple Hits')
                print(locus_tag)
                print(max_iden_dict)
                hit_list = ':'.join(max_iden_dict.keys())
                top_hits[locus_tag] = hit_list
    return top_hits


def identifyTopHits(blast_output_file):
    """
    :param blast_output_file: Output file from Blastp
    :return: Dictionary of top hits for each unannotated locus that passes the identity and coverage threshold
    """
    all_hits_dict = {}
    blast_output_raw = open(blast_output_file, 'r').readlines()
    for line in blast_output_raw:
        if line[0] == '#':
            continue
        line_elements = line.strip().split()
        if line_elements[0] not in all_hits_dict.keys():
            all_hits_dict[line_elements[0]] = {}
        iden = float(line_elements[5])
        if iden < args.identity:
            continue
        else:
            qcov = (float(line_elements[4])/float(line_elements[1])) * 100.0
            scov = (float(line_elements[4])/float(line_elements[3])) * 100.0
            if qcov >= args.coverage and scov >= args.coverage:
                corresponding_rv_hit = line_elements[2].split('|')[-1]
                all_hits_dict[line_elements[0]][corresponding_rv_hit] = (iden, qcov)
    top_hits_dict = filterTopHits(all_hits_dict)
    return top_hits_dict


def updateGenbankFile(gbk_file, top_hits_dict):
    """
    :param gbk_file: Genbank file that is to be updated
    :param top_hits_dict: Dictionary of top hits. The key is the 'L(2)_' locus tag and the value is corresponding top
    hit in H37Rv. If multiple top hits are present with same identity and coverage values, the value is all the top Rv
    hits separated by ':'
    :return: None. Updates the genbank files in place
    """
    genes_dict = genesLocusTagParser()
    genbank_filepath = args.annotation_directory + '/' + gbk_file
    modified_genbank = []
    for record in SeqIO.parse(genbank_filepath, 'genbank'):
        features_old = record.features[:]
        record.features = []
        for feature in features_old:
            if feature.type == 'CDS':
                locus_tag = feature.qualifiers['locus_tag'][0]
                if locus_tag in top_hits_dict.keys() and ':' not in top_hits_dict[locus_tag]:
                    # Locus tag has only one corresponding top hit
                    modified_locus_tag = top_hits_dict[locus_tag]
                    feature.qualifiers['locus_tag'] = [modified_locus_tag]
                    if modified_locus_tag in genes_dict.keys():
                        # If Rv has a corresponding gene name, this is updated in the feature as well
                        gene_name = genes_dict[modified_locus_tag]
                        feature.qualifiers['gene'] = [gene_name]
                elif locus_tag in top_hits_dict.keys() and ':' in top_hits_dict[locus_tag]:
                    # All corresponding hits are added ads a note
                    genes_list = top_hits_dict[locus_tag].split(':')
                    modified_locus_tag = genes_list[0]
                    feature.qualifiers['locus_tag'] = [modified_locus_tag]
                    if modified_locus_tag in genes_dict.keys():
                        gene_name = genes_dict[modified_locus_tag]
                        feature.qualifiers['gene'] = [gene_name]
                    all_hits = ','.join(genes_list[1:])
                    add_note = 'This gene also had positive hits to ' + all_hits
                    if 'note' in feature.qualifiers.keys():
                        feature.qualifiers['note'].append(add_note)
                    else:
                        feature.qualifiers['note'] = [add_note]
            record.features.append(feature)
        modified_genbank.append(record)
    gb_file_handle = open(genbank_filepath, 'w')
    SeqIO.write(modified_genbank, gb_file_handle, 'genbank')
    gb_file_handle.close()
    return


def main():

    # This script takes as options -i <isolate_id> -d <annotation_directory> -o <output_directory> -iden <identity> -cov
    # <coverage> from the commandline. The default values for the -i option is all, i.e. the script updates annotations
    # for all genbank files in the annotation_directory. The defaults for -d and -o options are
    # $GROUPHOME/data/annotation and None respectively. If -o option is not specified, intermediate output files with
    # top blast hits and distribution of amino acid % identity will not be printed/saved to file. The defaults for -iden
    # and -cov options are 85% and 99% respectively, where -id

    fp_annotation = os.environ['GROUPHOME'] + '/data/annotation'
    h37rv_protein_fasta = os.environ['GROUPHOME'] + '/resources/H37Rv-CDS-NC_000962.3.fasta'

    ###fp_annotation = 'test_update_annotation'
    ###h37rv_protein_fasta = 'test_update_annotation/H37Rv-CDS-NC_000962.3.fasta'
    purge_directory = True

    ####################################################################################################################
    ################################ Parsing Input options #############################################################
    ####################################################################################################################

    parser = argparse.ArgumentParser(description='Update annotation files from annomerge with H37Rv genes',
                                     epilog='Optional arguments --isolate, --annotation_directory, --output_directory, '
                                            '--identity and --coverage', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--isolate', help='Isolate ID (optional)', default='all')
    parser.add_argument('-d', '--annotation_directory', help='Path for annotation directory (optional)',
                        default=fp_annotation)
    parser.add_argument('-o', '--output_directory', help='If filepath for output_directory is specified, intermediate '
                                                         'output files with top blast hits and distribution of amino '
                                                         'acid percent identity will be printed/saved to path',
                        default=None)
    parser.add_argument('-iden', '--identity', help='Threshold for percentage amino acid sequence identity to a H37Rv '
                                                    'protein to be considered a confident hit. Annotation is '
                                                    'transferred from the top hit to H37Rv if sequence identity is '
                                                    'above this threshold. Default = 85', default=85)
    parser.add_argument('-cov', '--coverage', help='Threshold for sequence coverage of hit to H37Rv to be considered a '
                                                   'confident hit. Annotation is transferred from the top hit to H37Rv '
                                                   'if coverage of hit to H37Rv is above this threshold. Default = 99',
                        default=99)
    global args
    global output_dir_name
    args = parser.parse_args()
    try:
        directory_elements = os.listdir(args.annotation_directory)
        if len(directory_elements) == 0:
            sys.exit('Annotation directory is empty')
    except OSError:
        sys.exit('Annotation directory does not exist')
    files_in_anno_dir = os.listdir(args.annotation_directory)
    genbank_files = []
    if args.isolate == 'all':
        print('Isolate_id not specified. All isolates will be updated')
        for a_file in files_in_anno_dir:
            if a_file.endswith('.gbk') and a_file not in genbank_files:
                genbank_files.append(a_file)
    else:
        genbank_file = args.isolate + '.gbk'
        if genbank_file not in files_in_anno_dir:
            sys.exit('Genbank file does not exist for the given isolate')
        genbank_files.append(genbank_file)
    if args.output_directory is None:
        output_dir_name = args.annotation_directory + '/annotation_temp'
    else:
        output_dir_name = args.output_directory
    if os.path.isdir(output_dir_name) is True:
        purge_directory = False
        print('Warning: Output directory ' + output_dir_name + ' already exists. Temporary output files will be written'
                                                               ' to this directory and will not be purged.')
    else:
        os.mkdir(output_dir_name)
    try:
        identity_threshold = float(args.identity)
    except ValueError:
        sys.exit('Invalid amino acid sequence identity threshold')
    try:
        coverage_threshold = float(args.coverage)
    except ValueError:
        sys.exit('Invalid amino acid sequence coverage threshold')

    ####################################################################################################################
    ############### Getting unannotated CDSs, performing Blastp with H37Rv CDS and filtering positive hits #############
    ####################################################################################################################

    for gbk_file in genbank_files:
        isolate_id = gbk_file.split('.gbk')[0]
        print('Updating annotation for: ' + isolate_id)
        unannotated_genes_fasta = getUnannotatedGenes(gbk_file)
        blast_output = output_dir_name + '/unannotated_to_rv_blast'
        blast_command = NcbiblastpCommandline(query=unannotated_genes_fasta, subject=h37rv_protein_fasta,
                                              out=blast_output, outfmt='"7 qseqid qlen sseqid slen qlen length pident '
                                                                       'qcovs"')
        os.system(str(blast_command))
        top_hits = identifyTopHits(blast_output)
        updateGenbankFile(gbk_file, top_hits)
    if purge_directory:
        os.rmdir(output_dir_name)


if __name__ == "__main__":
    main()

