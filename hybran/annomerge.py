__author__ = "Deepika Gunasekaran"
__version__ = "2.0.0"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@sdsu.edu"
__status__ = "Development"

# Title: Merge annotation from Prokka for the positions which are not annotated by RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files in case of
# multiple contigs/chromosome annotations and a Genbank file (.gbk) file from Prokka annotation run with a reference and
# a Genbank file (.gbk) file from Prokka annotation run without a reference. The output is an EMBL file with annotation
# predominantly from RATT and the intergenic regions annotated by RATT are filled with Prokka. This script also
# generates a log file to indicate characteristics of the transferred features from Prokka.

import sys
import Bio
from Bio import SeqIO
from Bio.Seq import translate
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
import collections
from numpy import median
import os
import tempfile
import pickle
import logging
import time
import subprocess

from . import converter
from . import config


def load_reference_info(proteome_fasta):
    """
    This function gets all the global variables for annomerge.
    1. A list of genes from the reference protein fasta
    2. A list of locus tags from the reference protein fasta
    3. A dict with gene name as keys and locus tags as values
    4. A dict with locus tags as keys and gene name as values
    5. A dict with locus tags as keys and amino acid sequence as values
    6. A dict with locus tags as keys and path to temporary FASTA (with corresponding amino acid sequence) as values
    7. A dict with locus tags as keys and amino acid lengths as values
    8. A dict with locus tags as keys and nucleotide lengths as values
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    reference_gene_list = []
    reference_locus_list = []
    reference_locus_gene_dict = {}
    reference_gene_locus_dict = {}
    ref_temp_fasta_dict = {}
    ref_protein_lengths = {}

    ref_fasta_records = SeqIO.parse(proteome_fasta, 'fasta')
    for record in ref_fasta_records:
        gene_locus_name = record.description
        gene_locus = gene_locus_name.split(':')
        if len(gene_locus[0]) == 0:  # If locus tag is not specified
            gene = gene_locus[1]
            locus = gene
        elif len(gene_locus[1]) == 0:  # If gene name is not specified
            locus = gene_locus[0]
            gene = locus
        else:
            locus = gene_locus[0]
            gene = gene_locus[1]
        reference_gene_list.append(gene)
        reference_locus_list.append(locus)
        reference_gene_locus_dict[gene] = locus
        reference_locus_gene_dict[locus] = gene
        fp = tempfile.NamedTemporaryFile(suffix='.fasta',
                                         dir=hybran_tmp_dir,
                                         delete=False,
                                         mode='w')
        ref_temp_fasta_dict[locus] = fp.name
        header = '>' + locus + '\n'
        seq = str(record.seq)
        fp.write(header)
        fp.write(seq)
        ref_protein_lengths[locus] = len(seq)
    return reference_gene_list, reference_locus_list, reference_gene_locus_dict, reference_locus_gene_dict, \
           ref_temp_fasta_dict, ref_protein_lengths


def get_prom_for_gene(feature_list, source_seq):
    """
    This function gets the promoter sequence (40 bp upstream of gene) for a set of features, stores the
    sequences in a temporary nucleotide FASTA and generates a dictionary key to reference these FASTA by headers for
    each feature in the list
    :param feature_list: list of type SeqFeature (Biopython feature) formats
    :param source_seq: Nucleotide sequence of genome
    :return: Dictionary key to reference temporary FASTA file with promoter sequences, key is the locus tag from feature
     list and value is the FASTA header
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    ref_prom_dict = {}
    for f in feature_list:
        if f.type != 'CDS':
            continue
        feature_start = int(f.location.start)
        feature_stop = int(f.location.end)
        feature_strand = int(f.location.strand)
        locus = f.qualifiers['locus_tag'][0]
        if feature_strand == 1:
            if feature_start == 0:
                prom_end = len(source_seq)
                prom_start = len(source_seq) - 40
                prom_seq = str(source_seq[prom_start:prom_end])
            elif feature_start < 40:
                prom_end = feature_start
                prev_prom_len = 40 - prom_end
                prom_start = len(source_seq) - prev_prom_len
                prom_seq = str(source_seq[prom_start:len(source_seq)]) + str(source_seq[0:prom_end])
            else:
                prom_end = feature_start
                prom_start = feature_start - 40
                prom_seq = str(source_seq[prom_start:prom_end])
        else:
            prom_start = feature_stop
            prom_end = feature_stop + 40
            prom_seq = str(source_seq[prom_start:prom_end])
        fp_prom = tempfile.NamedTemporaryFile(suffix='_prom.fasta',
                                              dir=hybran_tmp_dir,
                                              delete=False,
                                              mode='w')
        ref_prom_dict[locus] = fp_prom.name
        header_prom = '>' + locus + '_prom\n'
        seq_prom = prom_seq + '\n'
        fp_prom.write(header_prom)
        fp_prom.write(seq_prom)
        fp_prom.close()
    return ref_prom_dict


def get_ordered_features(feature_list):
    """
    This function takes list of features and returns the list sorted by genomic location
    :param feature_list: list of type SeqFeature (Biopython feature) formats
    :return: sorted list of type SeqFeature (Biopython feature) formats, sorted by genomic location
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
    This function parses through a list of CDSs and returns a unique list of CDSs i.e. removes duplicate
    annotations which are in the same location in the genome with same locus tag
    :param feature_list: list of type SeqFeature (Biopython feature) formats
    :return: list of type SeqFeature (Biopython feature) formats after removing duplicate CDSs i.e. annotation of same
    gene in same position
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
    """
    This function takes as input a list of features and returns a dictionary with the key as a tuple of
    feature start and stop positions and the value as the feature.
    :param feature_list: List of features (SeqFeature objects)
    :return: sorted dictionary ordered by the genomic position i.e. feature location where key is a tuple
    (feature_start, feature_end, feature strand) and value is the corresponding SeqFeature object
    """
    feature_dict = {}
    for feature in feature_list:
        feature_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
        feature_dict[feature_key] = feature
    sorted_feature_dict = collections.OrderedDict(sorted(feature_dict.items()))
    return sorted_feature_dict


def check_for_dnaA(feature_list):
    """
    This function takes as input a list of features and checks if the genome was circularized and if dnaA
    is the first annotated CDS in the annotation file.
    :param feature_list: List of features in the RATT annotation record (SeqFeature objects)
    :return: None
    """

    logger = logging.getLogger('VerifyDnaA')
    logger.debug('Checking if dnaA is the first CDS')
    feature_dictionary = generate_feature_dictionary(feature_list)
    for cds in feature_dictionary.keys():
        if feature_dictionary[cds].type != 'CDS':
            continue
        # Checks for both 'Rv0001' if reference is H37Rv, else looks for dnaA in the gene_name
        if 'gene' not in feature_dictionary[cds].qualifiers.keys():
            gene_name = feature_dictionary[cds].qualifiers['locus_tag'][0]
        else:
            gene_name = feature_dictionary[cds].qualifiers['gene'][0]
        if feature_dictionary[cds].qualifiers['locus_tag'][0] != 'Rv0001' and gene_name != 'dnaA':
            logger.error('DnaA is not the first gene in this genome. Possible circularization error')
            logger.error('Exiting Annomerge. If you do not want dnaA check use the -illumina flag.')
            sys.exit()
        elif (feature_dictionary[cds].qualifiers['locus_tag'][0] == 'Rv0001' or gene_name == 'dnaA') and \
                int(feature_dictionary[cds].location.start) == 0:
            break
        else:
            logger.error('DnaA is the first gene in this genome. But the position is off.')
            logger.error('Exiting Annomerge. If you do not want dnaA check use the -illumina flag.')
            sys.exit()
    return


def get_ratt_corrected_genes(ratt_report_fp, reference_gene_list):
    """
    This function parses through the RATT results and gets the gene for which the start/stop coordinates
    were corrected by RATT
    :param ratt_report_fp: Correction Report from RATT
    :return: returns list of genes whose start/stop has been corrected by RATT
    """
    corrected_genes_list = []
    report_raw = open(ratt_report_fp).readlines()
    for line in report_raw:
        if len(line) <= 1:
            continue
        gene = line.strip().split()[0]
        if gene not in reference_gene_list:
            continue
        if gene not in corrected_genes_list:
            corrected_genes_list.append(gene)
    return corrected_genes_list


def rename_locus(gene, strand, reference_locus_list):
    """
    This function checks names of existing locus tags in the reference and names the newly merged gene
    ensuring that the new name does not conflict with existing locus tags
    :param gene: gene name to be modified
    :param strand: strand information to check if gene is present in complementary strand (- or '-1' or -1)
    :return: new gene name for merged gene
    """

    char_start = 65  # ASCII conversion of 'A'
    # CAUTION!!!
    # ASSUMTION: Format of the gene name/ locus tag from RATT
    gene_name = gene[:6]
    if strand == '-1' or strand == '-' or strand == -1:
        # Adding in 'c' to denote gene is located in the complementary strand
        new_gene_name = gene_name + chr(char_start) + 'c'
    else:
        new_gene_name = gene_name + chr(char_start)
    while new_gene_name in reference_locus_list:
        # If the assigned new gene name is in H37Rv, increment the alphabet suffix
        char_start += 1
        if strand == '-1' or strand == '-' or strand == -1:
            new_gene_name = gene_name + chr(char_start) + 'c'
        else:
            new_gene_name = gene_name + chr(char_start)
    return new_gene_name


def identify_merged_genes(ratt_features):
    """
    This function takes as input a list features annotated in RATT and identifies instances of merged CDS
    annotations. The function returns a Boolean value (True if such instances are identified) and a dictionary with
    the strand as the key and locations of the merged genes as values.
    :param ratt_features: List of features from RATT (SeqFeature objects)
    :return:
    Dictionary of merged genes with keys as strand (-1 or 1),
    list of gene location tuples (gene start, gene end) as values and boolean to indicate if merged genes are identified
        from RATT annotation.
    """
    logger = logging.getLogger('FindMergedGenes')
    ratt_annotations = {}   # This dictionary holds 2 keys '-1' and '+1' denoting the strand and locations of all CDS
    # annotations in RATT for each strand
    ratt_unmerged_genes = {}
    ratt_merged_genes = {-1: [], 1: []}
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
                logger.warning('UNACCOUNTED CASE')
                logger.warning(gene_location)
    return ratt_merged_genes, merged_genes


def get_annotation_for_merged_genes(merged_genes, prokka_features, ratt_features, reference_locus_list):
    """
    This function takes as input the dictionary of merged genes from the identify_merged_genes function,
    a list of Prokka features and a list of RATT features and returns annotation for the merged genes from,
    Prokka: If the genes are annotated as merged in both RATT and Prokka
    RATT: If Prokka does not annotate these genes as merged
    :param merged_genes: Dictionary of merged genes by strand (output from identify_merged_genes)
    :param prokka_features: List of features form Prokka (SeqFeature objects)
    :param ratt_features: List of features from RATT (SeqFeature objects)
    :return:
    """
    logger = logging.getLogger('AnnotateMergedGenes')
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
                                                  ' in reference) are merged in this isolate '
                                                  '(annotation from Prokka)']
                else:
                    feature.qualifiers['note'].append('The genes ' + merged_features_string +
                                                      ' in in reference) are merged in this isolate '
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
        new_gene_name = rename_locus(new_feature.qualifiers['locus_tag'][0],
                                     new_feature.location.strand,
                                     reference_locus_list)
        new_feature.qualifiers['locus_tag'] = [new_gene_name]
        if 'note' in new_feature.qualifiers.keys():
            new_feature.qualifiers['note'].append('The genes ' + merged_features_string +
                                                  ' in reference) are merged in this isolate '
                                                  '(annotation from RATT)')
        else:
            new_feature.qualifiers['note'] = ['The genes ' + merged_features_string +
                                              ' in reference) are merged in this isolate '
                                              '(annotation from RATT)']
        merged_features_addition.append(new_feature)
        merged_genes[new_feature.location.strand].remove(location)
    if len(merged_genes[-1]) > 0 or len(merged_genes[1]) > 0:
        corner_cases = True
    else:
        corner_cases = False
    return merged_features_addition, corner_cases, merged_genes, final_ratt_features, final_prokka_features


def blast_feature_sequence_to_ref(query, locus_tag, reference_locus_list, ref_temp_fasta_dict, seq_ident, seq_covg, to_print=False):
    """
    Blasts feature sequence to corresponding amino acid from reference FASTA file.
    :param query: Amino acid sequence (String format)
    :param locus_tag: Locus tag of corresponding amino acid sequence
    :param to_print: Boolean value to indicate whether BLAST results should be printed to log file
    :return: boolean value to indicate whether annotation should be added on not and dictionary of all BLAST hits with
    corresponding locus tag as key and list as values where list contains [percent identity of hit, subject coverage,
    query coverage]
    """

    logger = logging.getLogger('BlastFeatureToRef')
    if locus_tag in reference_locus_list:
        subject_fp = ref_temp_fasta_dict[locus_tag]
        blast_to_rv = NcbiblastpCommandline(subject=subject_fp, outfmt='"7 qseqid qlen sseqid slen qlen length'
                                                                       ' pident qcovs"')
        stdout, stderr = blast_to_rv(stdin=query)
        if to_print:
            logger.debug(stdout)
        hit, all_hits, note = identify_top_hits_mtb(stdout, identity=seq_ident, coverage=seq_covg)
        if hit == locus_tag:
            add_annotation = True
            hit_stats = all_hits[hit]
        else:
            add_annotation = False
            hit_stats = []
    else:
        add_annotation = True
        hit_stats = []
    return add_annotation, hit_stats


def isolate_valid_ratt_annotations(feature_list, ref_temp_fasta_dict, reference_locus_list, seq_ident, seq_covg):
    """
    This function takes as input a list of features and checks if the length of the CDSs are divisible by
    3 and if the CDS is split across multiple locations. If so, it outputs the features to stdout and removes them
    from the valid_ratt_annotations. The function BLASTs the sequence to corresponding amino acid sequence in
    reference as well.
    :param feature_list: List of features from RATT (list of SeqFeature objects)
    :return: List of valid RATT features (list of SeqFeature objects)
    """
    logger = logging.getLogger('ValidateRATTCDSs')
    logger.debug('Parsing through RATT annotations')
    unbroken_cds = []
    num_joins = 0
    non_cds_features = []
    broken_cds = []
    ratt_blast_results = {}
    for feature in feature_list:
        # Identify features with 'joins'
        if feature.type == 'CDS' and 'Bio.SeqFeature.CompoundLocation' in str(type(feature.location)):
            locus_tag = feature.qualifiers['locus_tag'][0]
            if locus_tag not in broken_cds:
                broken_cds.append(locus_tag)
            num_joins += 1
        elif feature.type == 'CDS' and feature.location is None:
            logger.warning('Invalid CDS: Location of CDS is missing')
            logger.warning(feature)
        elif feature.type == 'CDS' and (len(feature.location) % 3) != 0:
            logger.warning('Nucleotide sequence is not divisible by 3')
            logger.warning(feature)
        elif feature.type == 'CDS':
            unbroken_cds.append(feature)
        else:
            non_cds_features.append(feature)
    valid_features = []
    logger.debug("Valid CDSs before checking coverage: " + str(len(unbroken_cds)))
    for cds_feature in unbroken_cds:
        feature_sequence = translate(cds_feature.extract(record_sequence), table=11, to_stop=True)
        cds_locus_tag = cds_feature.qualifiers['locus_tag'][0]
        if len(feature_sequence) == 0:
            continue
        else:
            add_sequence, blast_stats = blast_feature_sequence_to_ref(str(feature_sequence), cds_locus_tag,
                                                                      reference_locus_list,
                                                                      ref_temp_fasta_dict,
                                                                      seq_ident=seq_ident,
                                                                      seq_covg=seq_covg)
            if add_sequence:
                if len(blast_stats) > 0:
                    ratt_blast_results[cds_locus_tag] = blast_stats
                cds_feature.qualifiers['translation'] = [str(feature_sequence)]
                valid_features.append(cds_feature)
            else:
                continue
    logger.debug("Valid CDSs after checking coverage: " + str(len(valid_features)))
    return valid_features, ratt_blast_results


def remove_duplicate_annotations(ratt_features, prokka_features_dictionary):
    """
    This function prunes and selects the features that are relevant in Prokka and discards features in
    Prokka that are annotated by RATT by taking into account the gene name and the position of the features
    :param ratt_features: List of features from RATT (list of SeqFeature objects)
    :param prokka_features_dictionary: sorted dictionary of features from Prokka ordered by the genomic position i.e.
    feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the corresponding
    SeqFeature object
    :return: sorted dictionary of features from Prokka that are not annotated by RATT ordered by the genomic position
    i.e. feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the
    corresponding SeqFeature object and dictionary of genes in RATT that overlap with prokka annotations.
    """

    prokka_features_not_in_ratt = prokka_features_dictionary.copy()
    ratt_overlapping_genes = {}
    for ratt_feature in ratt_features:
        if ratt_feature.type != 'CDS':
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
                    prokka_features_not_in_ratt.pop((prokka_start, prokka_end, prokka_strand), None)
                    prokka_duplicate_removed = True
                elif len(ratt_feature.location) == len(prokka_feature.location):
                    prokka_features_not_in_ratt.pop((prokka_start, prokka_end, prokka_strand), None)
                    prokka_duplicate_removed = True
        if not prokka_duplicate_removed and ratt_feature.type != 'mRNA':
            ratt_overlapping_genes[(ratt_start, ratt_end, ratt_strand)] = ratt_feature
    return prokka_features_not_in_ratt, ratt_overlapping_genes


def get_interregions(embl_record, intergene_length=1):
    """
    # Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
    # Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    # Do not remove this comment
    # This function was modified by Deepika Gunasekaran

    This function gets the genomic locations that do not have an coding-sequence (intergenic regions)
    :param embl_record: EMBL SeqRecord
    :param intergene_length: minimum length of integernic region (Default: 1)
    :return:
    SeqRecord of intergenic locations,
    list of intergenic positions where each element in the list is a tuple (start, end, strand),
    dictionary of genes preceding the intergenic regions where the key is a tuple (start of intergenic region, strand)
        and value is the SeqFeature preceding the intergenic region,
    dictionary of genes succeeding the intergenic regions where the key is a tuple (end of intergenic region, strand) and
        value is the SeqFeature succeeding the intergenic region
    """

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
            cds_list_minus.append((mystart, myend, -1))
            pre_intergene[(myend, -1)] = feature
            post_intergene[(mystart, -1)] = feature
        elif feature.strand == 1:
            cds_list_plus.append((mystart, myend, 1))
            pre_intergene[(myend, 1)] = feature
            post_intergene[(mystart, 1)] = feature
        else:
            sys.stderr.write("No strand indicated %d-%d. Assuming +\n" % (mystart, myend))
            cds_list_plus.append((mystart, myend, 1))
            pre_intergene[(myend, 1)] = feature
            post_intergene[(mystart, 1)] = feature
    cds_list_plus = sorted(cds_list_plus)
    cds_list_minus = sorted(cds_list_minus)
    for i, pospair in enumerate(cds_list_plus[1:]):
        # Compare current start position to previous end position
        last_end = cds_list_plus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end: this_start]
            strand_string = "+"
            intergenic_records.append(SeqRecord(intergene_seq,
                                                id="%s-ign-%d" % (seq_record.name, i),
                                                description="%s %d-%d %s" % (seq_record.name, last_end + 1,
                                                                             this_start, strand_string)))
            intergenic_positions.append((last_end + 1, this_start, strand_string))
    for i, pospair in enumerate(cds_list_minus[1:]):
        last_end = cds_list_minus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end: this_start]
            strand_string = "-"
            intergenic_records.append(SeqRecord(intergene_seq,
                                                id="%s-ign-%d" % (seq_record.name, i),
                                                description="%s %d-%d %s" % (seq_record.name, last_end + 1,
                                                                             this_start, strand_string)))
            intergenic_positions.append((last_end + 1, this_start, strand_string))
    return intergenic_records, intergenic_positions, pre_intergene, post_intergene


def check_inclusion_criteria(annotation_mapping_dict, embl_file, ratt_annotation, prokka_annotation, ratt_genes_check,
                             ratt_gene_location, ref_temp_fasta_dict, reference_gene_locus_dict,
                             reference_locus_gene_dict, ratt_blast_results, reference_locus_list, seq_ident, seq_covg):
    """
    This function compares RATT and Prokka annotations and resolves conflicting annotations

    :param annotation_mapping_dict:
    :param embl_file:
    :param ratt_annotation:
    :param prokka_annotation:
    :param ratt_genes_check:
    :param ratt_gene_location:
    :return:
    """

    logger = logging.getLogger('CheckInclusionCriteria')
    included = False
    # Check if feature types are the same. If not add feature to EMBL record
    if ratt_annotation.type != prokka_annotation.type:
        if prokka_annotation.type == 'CDS':
            mod_prokka_annotation, \
                invalid_prokka = validate_prokka_feature_annotation(prokka_annotation,
                                                                    prokka_noref_dictionary,
                                                                    ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                    reference_gene_locus_dict=reference_gene_locus_dict,
                                                                    reference_locus_gene_dict=reference_locus_gene_dict,
                                                                    ratt_blast_results=ratt_blast_results,
                                                                    seq_ident=seq_ident,
                                                                    seq_covg=seq_covg,
                                                                    check_ratt=True,
                                                                    ratt_features=ratt_genes_check,
                                                                    ratt_locations=ratt_gene_location,
                                                                    reference_locus_list=reference_locus_list)
        else:
            mod_prokka_annotation = prokka_annotation
        embl_file.features.append(mod_prokka_annotation)
        included = True
    # Check if gene names match and if they don't or if gene names are missing, keep both
    elif 'gene' in ratt_annotation.qualifiers.keys() and 'gene' in prokka_annotation.qualifiers.keys():
        if ratt_annotation.qualifiers['gene'] != prokka_annotation.qualifiers['gene']:
            mod_prokka_annotation, \
                invalid_prokka = validate_prokka_feature_annotation(prokka_annotation,
                                                                    prokka_noref_dictionary,
                                                                    ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                    reference_gene_locus_dict=reference_gene_locus_dict,
                                                                    reference_locus_gene_dict=reference_locus_gene_dict,
                                                                    ratt_blast_results=ratt_blast_results,
                                                                    seq_ident=seq_ident,
                                                                    seq_covg=seq_covg,
                                                                    check_ratt=True,
                                                                    ratt_features=ratt_genes_check,
                                                                    ratt_locations=ratt_gene_location,
                                                                    reference_locus_list=reference_locus_list)
            if not invalid_prokka:
                embl_file.features.append(mod_prokka_annotation)
                included = True
        # If gene names are the same and the lengths of the genes are comparable between RATT and Prokka annotation
        # (difference in length of less than/equal to 10 bps), the RATT annotation is preferred
        elif ratt_annotation.qualifiers['gene'] == prokka_annotation.qualifiers['gene'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 0:
            mod_prokka_annotation, \
                invalid_prokka = validate_prokka_feature_annotation(prokka_annotation,
                                                                    prokka_noref_dictionary,
                                                                    ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                    reference_gene_locus_dict=reference_gene_locus_dict,
                                                                    reference_locus_gene_dict=reference_locus_gene_dict,
                                                                    ratt_blast_results=ratt_blast_results,
                                                                    seq_ident=seq_ident,
                                                                    seq_covg=seq_covg,
                                                                    check_ratt=True,
                                                                    ratt_features=ratt_genes_check,
                                                                    ratt_locations=ratt_gene_location,
                                                                    reference_locus_list=reference_locus_list)
            if not invalid_prokka:
                embl_file.features.append(mod_prokka_annotation)
                included = True
    # If gene tag is missing and the product is not a hypothetical protein, check to see if the products are the
    # same between RATT and Prokka and if they are and if the lengths of the protein coding genes are comparable
    # preferred
    elif ('gene' not in ratt_annotation.qualifiers.keys() or 'gene' not in prokka_annotation.qualifiers.keys()) and \
            ('product' in ratt_annotation.qualifiers.keys() and 'product' in prokka_annotation.qualifiers.keys()):
        if ratt_annotation.qualifiers['product'] == 'hypothetical protein' or \
                        prokka_annotation.qualifiers['product'] == 'hypothetical protein':
            mod_prokka_annotation, \
                invalid_prokka = validate_prokka_feature_annotation(prokka_annotation,
                                                                    prokka_noref_dictionary,
                                                                    ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                    reference_gene_locus_dict=reference_gene_locus_dict,
                                                                    reference_locus_gene_dict=reference_locus_gene_dict,
                                                                    ratt_blast_results=ratt_blast_results,
                                                                    seq_ident=seq_ident,
                                                                    seq_covg=seq_covg,
                                                                    check_ratt=True,
                                                                    ratt_features=ratt_genes_check,
                                                                    ratt_locations=ratt_gene_location,
                                                                    reference_locus_list=reference_locus_list)
            if not invalid_prokka:
                embl_file.features.append(mod_prokka_annotation)
                included = True
        elif ratt_annotation.qualifiers['product'] == prokka_annotation.qualifiers['product'] and \
                        abs(len(prokka_annotation.location)-len(ratt_annotation.location)) > 0:
            mod_prokka_annotation, \
                invalid_prokka = validate_prokka_feature_annotation(prokka_annotation,
                                                                    prokka_noref_dictionary,
                                                                    ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                    reference_gene_locus_dict=reference_gene_locus_dict,
                                                                    reference_locus_gene_dict=reference_locus_gene_dict,
                                                                    ratt_blast_results=ratt_blast_results,
                                                                    seq_ident=seq_ident,
                                                                    seq_covg=seq_covg,
                                                                    check_ratt=True,
                                                                    ratt_features=ratt_genes_check,
                                                                    ratt_locations=ratt_gene_location,
                                                                    reference_locus_list=reference_locus_list)
            if not invalid_prokka:
                embl_file.features.append(mod_prokka_annotation)
                included = True
    else:
        logger.warning('CORNER CASE in check_inclusion_criteria')
    return embl_file, included


def get_top_hit(all_hits_dict):

    """
    This function parses through BLAST results to get the top hits.
    :param all_hits_dict: This dictionary consists of all Blast hits with corresponding unannotated locus_tag that have
    passed the defined thresholds for amino acid identity and coverage
    :return: returns only the elements of top hits with the locus tag. The key is the 'L(2)_' locus tag and the value is
     corresponding top hit in H37Rv. If multiple top hits are present with same identity and coverage values, the value
     is all the top Rv hits separated by ':'
    """

    top_hit_identity = 0
    for gene in all_hits_dict.keys():
        if all_hits_dict[gene][0] > top_hit_identity:
            top_hit_identity = all_hits_dict[gene]
            top_hit = gene
    return top_hit


def identify_top_hits_mtb(blast_output_file, identity=95, coverage=95, mtb=False):
    """
    This function gets the hits from BLAST tha passes the identity and coverage thresholds in the function
    :param blast_output_file: Output file from Blastp
    :param identity: identity threshold cutoff (Default: 95%)
    :param coverage: coverage threshold cutoff (Default: 95%)
    :param mtb: Internal parameter (to check for hits against previously reported MTB genes)
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
                all_hits_dict[corresponding_mtb_hit] = [iden, scov, qcov]
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
                all_hits_dict[corresponding_rv_hit] = [iden, scov, qcov]
            elif 85 <= iden < identity and qcov >= coverage and scov >= coverage:
                corresponding_rv_hit = line_elements[2].split('|')[-1]
                notes[corresponding_rv_hit] = [iden, scov, qcov]
            else:
                continue
        if rv_hit:
            top_hit = get_top_hit(all_hits_dict)
        elif len(list(all_hits_dict.keys())) != 0:
            top_hit = 'MTB:' + get_top_hit(all_hits_dict)
        else:
            top_hit = None
        return top_hit, all_hits_dict, notes


def get_essential_genes(gene_essentiality_fp=None):
    """
    This functions gets a set of Essential genes in H37Rv reported by DeJesus, 2017
    :param gene_essentiality_fp: Filepath to gene essentiality results (Default:
    'resources/gene_essentiality_dejesus_2017.tsv')
    :return: list od essential genes
    """
    # TODO : Checking for annotation validity against essential genes reported can be removed

    essential_genes = []
    if gene_essentiality_fp is None:
        gene_essentiality_fp = 'resources/gene_essentiality_dejesus_2017.tsv'
        gene_essentiality_raw = open(gene_essentiality_fp, 'r').readlines()
        for line in gene_essentiality_raw:
            line_elements = line.strip().split('\t')
            if line_elements[1] == 'ES':
                essential_genes.append(line_elements[0])
    else:
        gene_essentiality_raw = open(gene_essentiality_fp, 'r').readlines()
        essential_genes = [line.strip() for line in gene_essentiality_raw]
    return essential_genes


def validate_prokka_feature_annotation(feature, prokka_noref, reference_gene_locus_dict, reference_locus_gene_dict,
                                       ref_temp_fasta_dict, ratt_blast_results, reference_locus_list,
                                       seq_ident, seq_covg,
                                       check_ratt=False, ratt_features=[],
                                       ratt_locations={}):
    """

    :param feature:
    :param prokka_noref:
    :param check_ratt:
    :param ratt_features:
    :param ratt_locations:
    :param ratt_annotations:
    :return:
    """
    logger = logging.getLogger('ValidateProkkaFeatures')
    do_not_add_prokka = False
    mod_feature = feature
    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
    if feature.type != 'CDS':
        return mod_feature, do_not_add_prokka
    if 'gene' not in feature.qualifiers.keys():
        mod_feature = feature
    else:
        feature_seq = str(feature.qualifiers['translation'][0])
        if feature.qualifiers['gene'][0] in reference_locus_gene_dict.keys():
            locus_tag = feature.qualifiers['gene'][0]
            prokka_blast_list.append(locus_tag)
            retain_rv_annotation, blast_stats = blast_feature_sequence_to_ref(feature_seq, locus_tag,
                                                                              reference_locus_list=reference_locus_list,
                                                                              ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                              seq_ident=seq_ident,
                                                                              seq_covg=seq_covg)
            if retain_rv_annotation:
                mod_feature = feature
                if check_ratt:
                    if locus_tag in ratt_features and locus_tag not in ratt_blast_results.keys():
                        mod_feature = feature
                    elif locus_tag not in ratt_features:
                        mod_feature = feature
                    elif locus_tag in ratt_features and locus_tag in ratt_blast_results.keys():
                        ratt_start = int(list(ratt_locations.values())[0][0])
                        ratt_stop = int(list(ratt_locations.values())[0][1])
                        prom_mutation = False
                        if list(ratt_locations.values())[0][2] == 1:
                            if ratt_start == 0:
                                ratt_prom_end = len(record_sequence)
                                ratt_prom_start = len(record_sequence) - 40
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                            elif ratt_start < 40:
                                ratt_prom_end = ratt_start
                                prev_prom_len = 40 - ratt_prom_end
                                ratt_prom_start = len(record_sequence) - prev_prom_len
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:len(record_sequence)]) + \
                                                str(record_sequence[0:ratt_prom_end])
                            else:
                                ratt_prom_end = ratt_start
                                ratt_prom_start = ratt_start - 40
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                        else:
                            ratt_prom_start = ratt_stop
                            ratt_prom_end = ratt_stop + 40
                            ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                        blast_to_rv_prom = NcbiblastnCommandline(subject=ref_prom_fp_dict[locus_tag],
                                                                 outfmt='"7 qseqid sseqid pident length mismatch '
                                                                        'gapopen qstart qend sstart send evalue '
                                                                        'bitscore gaps"')
                        stdout, stderr = blast_to_rv_prom(stdin=str(ratt_prom_seq))
                        prom_blast_elements = stdout.split('\n')
                        for line in prom_blast_elements:
                            if line.startswith('#') or len(line) <= 1:
                                continue
                            blast_results = line.strip().split('\t')
                            if float(blast_results[2]) == 100.0 and int(blast_results[3]) == 40:
                                do_not_add_prokka = True
                            else:
                                prom_mutation = True
                        if prom_mutation is True:
                            ratt_coverage_measure = abs(int(ratt_blast_results[locus_tag][1]) -
                                                        int(ratt_blast_results[locus_tag][2]))
                            prokka_coverage_measure = abs(int(blast_stats[1]) - int(blast_stats[2]))
                            if ratt_coverage_measure < prokka_coverage_measure:
                                do_not_add_prokka = True
                            elif prokka_coverage_measure < ratt_coverage_measure:
                                logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                                if ratt_locations[locus_tag] in ratt_contig_features_dict.keys():
                                    ratt_contig_features_dict.pop(ratt_locations[locus_tag])
                                    popped_ratt_genes.append(locus_tag)
                                do_not_add_prokka = False
                            elif ratt_coverage_measure == prokka_coverage_measure:
                                # Checking for identity if coverage is the same
                                if int(ratt_blast_results[locus_tag][0]) > int(blast_stats[0]):
                                    do_not_add_prokka = True
                                elif int(blast_stats[0]) > int(ratt_blast_results[locus_tag][0]):
                                    logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                                    if ratt_locations[locus_tag] in ratt_contig_features_dict.keys():
                                        ratt_contig_features_dict.pop(ratt_locations[locus_tag])
                                        popped_ratt_genes.append(locus_tag)
                                    do_not_add_prokka = False
                                else:
                                    # If RATT and Prokka annotations are same, choose RATT
                                    do_not_add_prokka = True
                            else:
                                do_not_add_prokka = False
            else:
                if loc_key in prokka_noref.keys():
                    mod_feature = prokka_noref[loc_key]
                    if 'gene' in mod_feature.qualifiers.keys():
                        if 'gene_synonym' in mod_feature.qualifiers.keys():
                            mod_feature.qualifiers['gene_synonym'].append(mod_feature.qualifiers['gene'])
                        else:
                            mod_feature.qualifiers['gene_synonym'] = mod_feature.qualifiers['gene']
                    mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
                else:
                    mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
        elif feature.qualifiers['gene'][0] in reference_gene_locus_dict.keys():
            locus_tag = reference_gene_locus_dict[feature.qualifiers['gene'][0]]
            prokka_blast_list.append(locus_tag)
            retain_rv_annotation, blast_stats = blast_feature_sequence_to_ref(feature_seq, locus_tag,
                                                                              reference_locus_list=reference_locus_list,
                                                                              ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                              seq_ident=seq_ident,
                                                                              seq_covg=seq_covg)
            if retain_rv_annotation:
                mod_feature = feature
                if check_ratt:
                    if locus_tag in ratt_features and locus_tag not in ratt_blast_results.keys():
                        mod_feature = feature
                    elif locus_tag not in ratt_features:
                        mod_feature = feature
                    elif locus_tag in ratt_features and locus_tag in ratt_blast_results.keys():
                        ratt_start = int(list(ratt_locations.values())[0][0])
                        ratt_stop = int(list(ratt_locations.values())[0][1])
                        prom_mutation = False
                        if list(ratt_locations.values())[0][2] == 1:
                            if ratt_start == 0:
                                ratt_prom_end = len(record_sequence)
                                ratt_prom_start = len(record_sequence) - 40
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                            elif ratt_start < 40:
                                ratt_prom_end = ratt_start
                                prev_prom_len = 40 - ratt_prom_end
                                ratt_prom_start = len(record_sequence) - prev_prom_len
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:len(record_sequence)]) + \
                                                str(record_sequence[0:ratt_prom_end])
                            else:
                                ratt_prom_end = ratt_start
                                ratt_prom_start = ratt_start - 40
                                ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                        else:
                            ratt_prom_start = ratt_stop
                            ratt_prom_end = ratt_stop + 40
                            ratt_prom_seq = str(record_sequence[ratt_prom_start:ratt_prom_end])
                        blast_to_rv_prom = NcbiblastnCommandline(subject=ref_prom_fp_dict[locus_tag],
                                                                 outfmt='"7 qseqid sseqid pident length mismatch '
                                                                        'gapopen qstart qend sstart send evalue '
                                                                        'bitscore gaps"')
                        stdout, stderr = blast_to_rv_prom(stdin=str(ratt_prom_seq))
                        prom_blast_elements = stdout.split('\n')
                        for line in prom_blast_elements:
                            if line.startswith('#') or len(line) <= 1:
                                continue
                            blast_results = line.strip().split('\t')
                            if float(blast_results[2]) == 100.0 and int(blast_results[3]) == 40:
                                do_not_add_prokka = True
                            else:
                                prom_mutation = True
                        if prom_mutation is True:
                            ratt_coverage_measure = abs(int(ratt_blast_results[locus_tag][1]) -
                                                        int(ratt_blast_results[locus_tag][2]))
                            prokka_coverage_measure = abs(int(blast_stats[1]) - int(blast_stats[2]))
                            if ratt_coverage_measure < prokka_coverage_measure:
                                do_not_add_prokka = True
                            elif prokka_coverage_measure < ratt_coverage_measure:
                                logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                                if ratt_locations[locus_tag] in ratt_contig_features_dict.keys():
                                    ratt_contig_features_dict.pop(ratt_locations[locus_tag])
                                    popped_ratt_genes.append(locus_tag)
                                do_not_add_prokka = False
                            elif ratt_coverage_measure == prokka_coverage_measure:
                                # Checking for identity if coverage is the same
                                if int(ratt_blast_results[locus_tag][0]) > int(blast_stats[0]):
                                    do_not_add_prokka = True
                                elif int(blast_stats[0]) > int(ratt_blast_results[locus_tag][0]):
                                    logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                                    if ratt_locations[locus_tag] in ratt_contig_features_dict.keys():
                                        ratt_contig_features_dict.pop(ratt_locations[locus_tag])
                                        popped_ratt_genes.append(locus_tag)
                                    do_not_add_prokka = False
                                else:
                                    # If RATT and Prokka annotations are same, choose RATT
                                    do_not_add_prokka = True
                            else:
                                do_not_add_prokka = False
            else:
                if loc_key in prokka_noref.keys():
                    mod_feature = prokka_noref[loc_key]
                    if 'gene' in mod_feature.qualifiers.keys():
                        if 'gene_synonym' in mod_feature.qualifiers.keys():
                            mod_feature.qualifiers['gene_synonym'].append(mod_feature.qualifiers['gene'])
                        else:
                            mod_feature.qualifiers['gene_synonym'] = mod_feature.qualifiers['gene']
                    mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
                else:
                    mod_feature.qualifiers['gene'] = mod_feature.qualifiers['locus_tag']
        else:
            mod_feature = feature
    return mod_feature, do_not_add_prokka


def blast_ratt_anno_with_prokka_ltag(ratt_overlap_feature, prokka_overlap_feature, seq_ident, seq_covg):
    """

    :param ratt_overlap_feature:
    :param prokka_overlap_feature:
    :return:
    """

    hybran_tmp_dir = config.hybran_tmp_dir
    same_gene = False
    fp = tempfile.NamedTemporaryFile(suffix='_ratt.fasta',
                                     dir=hybran_tmp_dir,
                                     delete=False,
                                     mode='w')
    header = '>' + str(ratt_overlap_feature.qualifiers['locus_tag'][0]) + '\n'
    seq = str(ratt_overlap_feature.qualifiers['translation'][0])
    fp.write(header)
    fp.write(seq)
    fp.close()
    prokka_seq = str(prokka_overlap_feature.qualifiers['translation'][0])
    try:
        blast_to_ratt_rv = NcbiblastpCommandline(subject=fp.name, outfmt='"7 qseqid qlen sseqid slen qlen length pident'
                                                                         ' qcovs"')
        stdout_rv, stderr_rv = blast_to_ratt_rv(stdin=prokka_seq)
        is_hit = True
    except Bio.Application.ApplicationError:
        is_hit = False
    if is_hit:
        hit, all_hits, note = identify_top_hits_mtb(stdout_rv, identity=seq_ident, coverage=seq_covg)
        if len(all_hits) == 0:
            same_gene = False
        else:
            same_gene = True
    return same_gene


def correct_start_coords_prokka(prokka_record, correction_dict, fasta_seq, rv_seq, rv_cds_dict, reference_locus_list,
                                reference_gene_locus_dict):
    """
    This function parses through prokka records and corrects start coordinates for cases where Prodigal
    annotates these incorrectly
    :param prokka_record:
    :param correction_dict:
    :param fasta_seq:
    :param rv_seq:
    :param rv_cds_dict:
    :return:
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('CorrectCoords')
    modified_prokka_record = prokka_record[:]
    modified_prokka_record.annotations['comment'] = 'Merged reference based annotation from RATT and ab initio ' \
                                                    'annotation from Prokka'
    modified_prokka_record.features = []
    modified_features = []
    # Write Rv genes to temp file
    rv_temp_nuc_fasta_dict = {}
    rv_temp_prom_fasta_dict = {}
    for rv in correction_dict.keys():
        rv_feature = rv_cds_dict[rv]
        rv_feat_seq = str(rv_seq)[int(rv_feature.location.start):int(rv_feature.location.end)]
        rv_prom_seq = str(rv_seq)[int(rv_feature.location.start)-40:int(rv_feature.location.start)]
        fp = tempfile.NamedTemporaryFile(suffix='_nuc.fasta',
                                         dir=hybran_tmp_dir,
                                         delete=False,
                                         mode='w')
        fp_prom = tempfile.NamedTemporaryFile(suffix='_prom.fasta',
                                              dir=hybran_tmp_dir,
                                              delete=False,
                                              mode='w')
        rv_temp_nuc_fasta_dict[rv] = fp.name
        rv_temp_prom_fasta_dict[rv] = fp_prom.name
        header = '>' + rv + '\n'
        seq = rv_feat_seq + '\n'
        fp.write(header)
        fp.write(seq)
        header_prom = '>' + rv + '_prom\n'
        seq_prom = rv_prom_seq + '\n'
        fp_prom.write(header_prom)
        fp_prom.write(seq_prom)
        fp.close()
        fp_prom.close()
    for feature_prokka in prokka_record.features:
        if feature_prokka.type != 'CDS':
            modified_features.append(feature_prokka)
            continue
        if feature_prokka.type == 'CDS' and 'gene' not in feature_prokka.qualifiers.keys():
            rv_id = ''
        elif feature_prokka.type == 'CDS' and 'gene' in feature_prokka.qualifiers.keys():
            if feature_prokka.qualifiers['gene'][0] in reference_locus_list:
                rv_id = feature_prokka.qualifiers['gene'][0]
            else:
                if '_' in feature_prokka.qualifiers['gene'][0]:
                    gene_name = feature_prokka.qualifiers['gene'][0].split('_')[0]
                else:
                    gene_name = feature_prokka.qualifiers['gene'][0]
                if gene_name in reference_gene_locus_dict.keys():
                    rv_id = reference_gene_locus_dict[gene_name]
                else:
                    rv_id = ''
        else:
            rv_id = ''
        if feature_prokka.type == 'CDS' and len(rv_id) > 0:
            original_start = int(feature_prokka.location.start)
            original_end = int(feature_prokka.location.end)
            original_strand = int(feature_prokka.location.strand)
            original_location = FeatureLocation(ExactPosition(original_start), ExactPosition(original_end),
                                                strand=original_strand)
            if rv_id in correction_dict.keys():
                # If Prokka annotation of gene is in a different strand, prodigal start coordinates are considered
                # as is
                query_temp = tempfile.NamedTemporaryFile(suffix='_query_nuc.fasta',
                                                         dir=hybran_tmp_dir,
                                                         delete=False,
                                                         mode='w')
                prokka_nuc_seq = str(fasta_seq)[int(feature_prokka.location.start):int(feature_prokka.location.end)]
                q_header = '>' + rv_id + '_query\n'
                query_temp.write(q_header)
                query_temp.write(prokka_nuc_seq)
                query_temp.write('\n')
                blast_to_rv = NcbiblastnCommandline(subject=rv_temp_nuc_fasta_dict[rv_id],
                                                    outfmt='"7 qseqid sseqid pident length mismatch gapopen qstart qend'
                                                           ' sstart send evalue bitscore"')
                stdout, stderr = blast_to_rv(prokka_nuc_seq)
                stdout_elements = stdout.split('\n')

                blast_length = len(stdout_elements)
                blast_counter = 0
                for line in stdout_elements:
                    if line.startswith('#') or len(line) == 0:
                        blast_counter += 1
                        continue
                    else:
                        line_elements = line.strip().split('\t')
                        if int(line_elements[6]) < int(line_elements[8]):
                            if int(feature_prokka.location.start) < 40 or int(rv_cds_dict[rv_id].location.start) < 40:
                                change_start = int(line_elements[8]) - int(line_elements[6])
                                mod_feature = feature_prokka
                                mod_start = int(feature_prokka.location.start) - change_start
                                mod_end = int(feature_prokka.location.end)
                                mod_strand = int(feature_prokka.location.strand)
                                mod_feature.location = FeatureLocation(ExactPosition(mod_start), ExactPosition(mod_end),
                                                                       strand=mod_strand)
                                mod_feature_seq = translate(mod_feature.extract(fasta_seq), table=11, to_stop=True)
                                check_feature_seq = translate(mod_feature.extract(fasta_seq), table=11, to_stop=False)
                                if len(mod_feature_seq) == len(check_feature_seq) - 1:
                                    mod_feature.qualifiers['translation'] = [str(mod_feature_seq)]
                                    modified_features.append(mod_feature)
                                else:
                                    feature_prokka.location = original_location
                                    modified_features.append(feature_prokka)
                            else:
                                change_start = int(line_elements[8]) - int(line_elements[6])
                                rv_feature = rv_cds_dict[rv_id]
                                rv_prom_end = int(rv_feature.location.start)
                                rv_prom_start = int(rv_feature.location.start) - 40
                                rv_prom_location = FeatureLocation(ExactPosition(rv_prom_start),
                                                                   ExactPosition(rv_prom_end),
                                                                   strand=int(rv_feature.location.strand))
                                rv_prom_seq = rv_prom_location.extract(rv_seq)
                                mod_prokka_start = int(feature_prokka.location.start) - change_start
                                prokka_prom_start = mod_prokka_start - 40
                                prokka_prom_location = FeatureLocation(ExactPosition(prokka_prom_start),
                                                                       ExactPosition(mod_prokka_start),
                                                                       strand=int(feature_prokka.location.strand))
                                feature_prokka_prom_seq = prokka_prom_location.extract(fasta_seq)
                                blast_to_rv_prom = NcbiblastnCommandline(subject=rv_temp_prom_fasta_dict[rv_id],
                                                                         outfmt='"7 qseqid sseqid pident length '
                                                                                'mismatch gapopen qstart qend sstart '
                                                                                'send evalue bitscore gaps"')
                                stdout_2, stderr_2 = blast_to_rv_prom(stdin=str(feature_prokka_prom_seq))
                                prom_blast_elements = stdout_2.split('\n')
                                prom_blast = len(prom_blast_elements)
                                prom_counter = 0
                                for line in prom_blast_elements:
                                    if line.startswith('#') or len(line) == 0:
                                        prom_counter += 1
                                        continue
                                    else:
                                        line_elements = line.strip().split('\t')
                                        if float(line_elements[2]) == 100.0 and int(line_elements[12]) == 0:
                                            mod_feature = feature_prokka
                                            mod_start = int(feature_prokka.location.start) - change_start
                                            mod_end = int(feature_prokka.location.end)
                                            mod_strand = int(feature_prokka.location.strand)
                                            mod_feature.location = FeatureLocation(ExactPosition(mod_start), ExactPosition(mod_end), strand=mod_strand)
                                            mod_feature_seq = translate(mod_feature.extract(fasta_seq), table=11,
                                                                        to_stop=True)
                                            check_feature_seq = translate(mod_feature.extract(fasta_seq), table=11,
                                                                          to_stop=False)
                                            if len(mod_feature_seq) == len(check_feature_seq) - 1:
                                                mod_feature.qualifiers['translation'] = [str(mod_feature_seq)]
                                                modified_features.append(mod_feature)
                                            else:
                                                feature_prokka.location = original_location
                                                modified_features.append(feature_prokka)
                                        else:
                                            modified_features.append(feature_prokka)
                                if prom_blast == prom_counter:
                                    modified_features.append(feature_prokka)
                        elif int(line_elements[6]) > int(line_elements[8]):
                            if int(feature_prokka.location.start) < 40 or int(rv_cds_dict[rv_id].location.start) < 40:
                                change_start = int(line_elements[6]) - int(line_elements[8])
                                mod_feature = feature_prokka
                                mod_start = int(feature_prokka.location.start) + change_start
                                mod_end = int(feature_prokka.location.end)
                                mod_strand = int(feature_prokka.location.strand)
                                mod_feature.location = FeatureLocation(ExactPosition(mod_start),
                                                                       ExactPosition(mod_end), strand=mod_strand)
                                mod_feature_seq = translate(mod_feature.extract(fasta_seq), table=11, to_stop=True)
                                check_feature_seq = translate(mod_feature.extract(fasta_seq), table=11, to_stop=False)
                                if len(mod_feature_seq) == len(check_feature_seq) - 1:
                                    mod_feature.qualifiers['translation'] = [str(mod_feature_seq)]
                                    modified_features.append(mod_feature)
                                else:
                                    feature_prokka.location = original_location
                                    modified_features.append(feature_prokka)
                            else:
                                change_start = int(line_elements[6]) - int(line_elements[8])
                                rv_feature = rv_cds_dict[rv_id]
                                rv_prom_end = int(rv_feature.location.start)
                                rv_prom_start = int(rv_feature.location.start) - 40
                                rv_prom_location = FeatureLocation(ExactPosition(rv_prom_start),
                                                                   ExactPosition(rv_prom_end),
                                                                   strand=int(rv_feature.location.strand))
                                rv_prom_seq = rv_prom_location.extract(rv_seq)
                                mod_prokka_start = int(feature_prokka.location.start) + change_start
                                prokka_prom_start = mod_prokka_start - 40
                                prokka_prom_location = FeatureLocation(ExactPosition(prokka_prom_start),
                                                                       ExactPosition(mod_prokka_start),
                                                                       strand=int(feature_prokka.location.strand))
                                feature_prokka_prom_seq = prokka_prom_location.extract(fasta_seq)
                                blast_to_rv_prom = NcbiblastnCommandline(subject=rv_temp_prom_fasta_dict[rv_id],
                                                                         outfmt='"7 qseqid sseqid pident length '
                                                                                'mismatch gapopen qstart qend sstart '
                                                                                'send evalue bitscore gaps"')
                                stdout_2, stderr_2 = blast_to_rv_prom(stdin=str(feature_prokka_prom_seq))
                                prom_blast_elements = stdout_2.split('\n')
                                prom_blast = len(prom_blast_elements)
                                prom_counter = 0
                                for line in prom_blast_elements:
                                    if line.startswith('#') or len(line) == 0:
                                        prom_counter += 1
                                        continue
                                    else:
                                        line_elements = line.strip().split('\t')
                                        if float(line_elements[2]) == 100.0 and int(line_elements[12]) == 0:
                                            mod_feature = feature_prokka
                                            mod_start = int(feature_prokka.location.start) + change_start
                                            mod_end = int(feature_prokka.location.end)
                                            mod_strand = int(feature_prokka.location.strand)
                                            mod_feature.location = FeatureLocation(ExactPosition(mod_start),
                                                                                   ExactPosition(mod_end),
                                                                                   strand=mod_strand)
                                            mod_feature_seq = translate(mod_feature.extract(fasta_seq), table=11,
                                                                        to_stop=True)
                                            check_feature_seq = translate(mod_feature.extract(fasta_seq), table=11,
                                                                          to_stop=False)
                                            if len(mod_feature_seq) == len(check_feature_seq) - 1:
                                                mod_feature.qualifiers['translation'] = [str(mod_feature_seq)]
                                                modified_features.append(mod_feature)
                                            else:
                                                feature_prokka.location = original_location
                                                modified_features.append(feature_prokka)
                                        else:
                                            modified_features.append(feature_prokka)
                                if prom_blast == prom_counter:
                                    modified_features.append(feature_prokka)
                        else:
                            # The start and stop are the same in both Blast and H37Rv i.e. Prodigal start location is
                            # correct
                            modified_features.append(feature_prokka)
                if blast_length == blast_counter:
                    modified_features.append(feature_prokka)
            else:
                modified_features.append(feature_prokka)
        else:
            modified_features.append(feature_prokka)
    modified_prokka_record.features = get_ordered_features(modified_features)
    return modified_prokka_record


def get_feature_by_locustag(rec):
    """

    :param rec:
    :return:
    """
    global ref_genes_positions
    ref_genes_positions = {}
    cds_dict = {}
    for feature in rec.features:
        if feature.type != 'CDS':
            continue
        else:
            rv = feature.qualifiers['locus_tag'][0]
            cds_dict[rv] = feature
            ref_genes_positions[rv] = (int(feature.location.start), int(feature.location.end))
    return cds_dict


def blast(seq1, seq2):
    """

    AUTHOR: Sarah Ramirez-Busby
    :param seq1:
    :param seq2:
    :return:
    """
    blast_to_rv = NcbiblastnCommandline(subject=seq1,
                                        outfmt='"7 qseqid sseqid pident length mismatch '
                                               'gapopen qstart qend sstart send evalue bitscore"')
    stdout, stderr = blast_to_rv(str(seq2))
    stdout_elements = stdout.split('\n')
    return stdout_elements


def pick_best_hit(ratt_feature, prokka_feature, isolate_sequence):
    """

    AUTHOR: Sarah Ramirez-Busby
    MODIFIED BY: Deepika Gunasekaran
    :param ratt_feature:
    :param prokka_feature:
    :return:
    """

    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('BestFeatureHit')
    gene = ratt_feature.qualifiers['locus_tag'][0]
    if gene not in ref_genes_positions.keys():
        take_ratt = True
        return take_ratt
    startpos = ref_genes_positions[gene][0]
    stoppos = ref_genes_positions[gene][1]
    prokka_seq = isolate_sequence[prokka_feature.location.start:prokka_feature.location.end]
    ratt_seq = isolate_sequence[ratt_feature.location.start:ratt_feature.location.end]
    ref_gene_seq = ref_sequence[startpos - 1: stoppos]

    f = tempfile.NamedTemporaryFile(prefix=gene,
                                     suffix='.fasta',
                                     dir=hybran_tmp_dir,
                                     delete=False, mode='w')
    gene_fasta = f.name
    SeqIO.write(SeqRecord(ref_gene_seq, id=gene, description=''), gene_fasta, 'fasta')

    blast_out = blast(gene_fasta, ratt_seq)
    blast_results = False
    for i in blast_out:
        if i.startswith('Query'):
            sstart = int(i.split('\t')[8])
            send = int(i.split('\t')[9])
            blast_results = True
    if blast_results:
        # If the gene is positive strand
        if ratt_feature.location.strand == '+':
            # The start of the prokka gene minus where the prokka sequence aligns to H37Rv plus one to include the
            # position
            ratt_prom_end = ratt_feature.location.start.position - sstart + 1
            ratt_prom_start = ratt_prom_end - 40
            rv_prom_end = startpos - 40
            rv_seq = ref_sequence[startpos - 1: rv_prom_end]
        # If the gene is negative strand
        else:
            rv_prom_end = stoppos + 40
            # Where the prokka sequence aligns to the H37Rv seq plus the end coordinate of the prokka gene
            ratt_prom_start = send + ratt_feature.location.end.position
            ratt_prom_end = ratt_prom_start + 40
            rv_seq = ref_sequence[stoppos - 1: rv_prom_end]

        prom_f = tempfile.NamedTemporaryFile(prefix=gene,
                                         suffix='-prom.fasta',
                                         dir=hybran_tmp_dir,
                                         delete=False, mode='w')
        prom_fasta = prom_f.name
        SeqIO.write(SeqRecord(rv_seq, id=gene, description=''), prom_fasta, 'fasta')

        take_ratt = True
        prom_blast = blast(prom_fasta,
                           isolate_sequence[ratt_prom_start:ratt_prom_end])
        for i in prom_blast:
            if i.startswith('Query'):
                if int(i.split('\t')[2]) < 100.0:
                    take_ratt = False
                else:
                    take_ratt = True
    else:
        take_ratt = False
    return take_ratt


def fix_embl_id_line(embl_file):
    """

    :param embl_file:
    :return:
    """
    lines = []
    with open(embl_file, 'r') as embl:
        for line in embl:
            if line.startswith('ID'):
                lines.append(line.replace(' ; ; ; ; ; ', ' ; ; ; ; ; ; '))
            else:
                lines.append(line)
    with open(embl_file, 'w') as out:
        for line in lines:
            out.write(line)


def run_prodigal(reference_genome, temp_dir):
    """
    :param reference_genome:
    :param temp_dir: str path to Hybran temporary directory
    :return:
    """
    logger = logging.getLogger('Prodigal')
    c = os.getcwd()
    logger.debug('Executing Prodigal on ' + reference_genome)
    cmd = [os.sep.join([script_dir, 'prodigal.sh']),
           reference_genome,
           temp_dir]
    subprocess.call(cmd)
    os.chdir(c)


def run(isolate_id, annotation_fp, ref_proteins_fasta, ref_embl_fp, reference_genome, script_directory, seq_ident, seq_covg, illumina=False,
        fill_gaps=True, check_mtb=False, mtb_fasta_fp='', essentiality=False):
    """
    Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file> -m
    <output_merged_genes> from the commandline. The log file output stats about the features that are added to the
    RATT annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/'isolate_id'.gbk and
    'isolate_id'/annomerge/'isolate_id'.log

    :param isolate_id: ID of the isolate (Example: H37Rv, 1-0006, etc.). This is the isolate_id that is used for naming
     Genbank files in Prokka
    :param annotation_fp: Filepath where RATT, Prokka, reference and prokka no-reference annotations are located.
    Annomerge assumes that RATT annotations are located in <annotation_fp>/ratt, Prokka reference annotations are
    located in <annotation_fp>/prokka and prokka annotations without reference is located in
    <annotation_fp>/prokka-noreference. Additionally annomerge also assumes that withing prokka and prokka-noreference
    directories, the genbank files are located in <isolate_id>.gbk
    :param ref_proteins_fasta: File path for proteome fasta of reference strain
    :param ref_embl_fp: File path for annotated EMBL file for reference strain
    :param reference_genome: File path for nucleotide fasta of assembled genome
    :param script_dir: Directory where hybran scripts are located
    :param illumina: Flag to check circularization. Default is False. If set to True, circularization with not be
    checked and dnaA might not be the first gene.
    :param fill_gaps: Flag to fill gaps in annotation from prokka no-reference. Default is true. If set to false,
    unannotated regions will not be check against prokka noreference run.
    :param check_mtb: Flag to check if genes should be blasted against existing set of novel genes. Default is false
    and if set to True, mtb_fasta_fp should be provided.
    :param mtb_fasta_fp: File path for amino acid fasta of novel genes.
    :param essentiality: Check annotation of essential genes in isolate
    :return: EMBL record (SeqRecord) of annotated isolate
    """

    hybran_tmp_dir = config.hybran_tmp_dir
    global script_dir
    script_dir = script_directory
    logger = logging.getLogger('Annomerge')
    logger.debug('Running Annomerge on ' + isolate_id)
    start_time = time.time()

    reference_gene_list, reference_locus_list, reference_gene_locus_dict, reference_locus_gene_dict, \
        ref_temp_fasta_dict, ref_protein_lengths = load_reference_info(ref_proteins_fasta)
    if annotation_fp.endswith('/'):
        file_path = annotation_fp + isolate_id + '/'
    else:
        file_path = annotation_fp + '/' + isolate_id + '/'
    ratt_file_path = file_path + 'ratt'
    ratt_correction_files = []
    ratt_gbk_files = []
    try:
        ratt_embl_files = [embl_file for embl_file in os.listdir(ratt_file_path) if embl_file.endswith('.final.embl')]
        for embl_file in ratt_embl_files:
            gbk = converter.convert_embl_to_gbk(ratt_file_path + '/' + embl_file)
            ratt_gbk_files.append(gbk)
        correction_files = [cf for cf in os.listdir(ratt_file_path) if cf.endswith('.Report.txt')]
        for corr_file in correction_files:
            corr_file_path = ratt_file_path + '/' + corr_file
            ratt_correction_files.append(corr_file_path)
    except OSError:
        logger.error('Expecting RATT annotation files but found none')
    if not ratt_gbk_files:
        logger.error('RATT did not complete running. Please see the log for more details.')
    try:
        input_prokka_genbank = file_path + 'prokka/' + isolate_id + '.gbk'
    except OSError:
        logger.error('Expecting Prokka annotation file but found none')
    output_merged_genes = isolate_id + '/annomerge/merged_genes.gbk'
    output_genbank = isolate_id + '.gbk'
    add_noref_annotations = fill_gaps
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    isolate_sequence = prokka_records[0].seq
    prokka_record_fp = file_path + 'prokka-noreference/' + isolate_id + '.gbk'
    prokka_record_noref = list(SeqIO.parse(prokka_record_fp, 'genbank'))
    annomerge_records = []

    # RUNNING PRODIGAL
    embl_dict = {}
    prodigal_list = []
    incorrect_coords_dict = {}
    prodigal_results_fp = hybran_tmp_dir + '/prodigal-test/reference_prodigal'
    if not os.path.exists(prodigal_results_fp):
        run_prodigal(reference_genome, hybran_tmp_dir)

    dict_save_fp = prodigal_results_fp + 'reference_prodigal.p'
    ref_embl_record = SeqIO.read(ref_embl_fp, 'embl')
    prodigal_raw = open(prodigal_results_fp, "r").readlines()
    for feat in ref_embl_record.features:
        if feat.type != 'CDS':
            continue
        embl_dict[(int(feat.location.start) + 1, int(feat.location.end))] = feat.qualifiers['locus_tag'][0]
    ordered_embl_dict = collections.OrderedDict(sorted(embl_dict.items()))
    for line in prodigal_raw:
        if line.startswith('     CDS'):
            line_elements = line.strip().split()
            coords = line_elements[1]
            if '(' in coords:
                coords_elements = coords.split('..')
                coords_start = coords_elements[0].split('(')[1]
                coords_end = coords_elements[1].split(')')[0]
                prodigal_list.append((int(coords_start), int(coords_end)))
            else:
                coords_elements = coords.split('..')
                prodigal_list.append((int(coords_elements[0].replace('<', '')), int(coords_elements[1])))
    for gene_coord in prodigal_list:
        for gene in ordered_embl_dict.keys():
            if gene[0] > gene_coord[1]:
                break
            elif gene[1] == gene_coord[1] and gene[0] != gene_coord[0]:
                prodigal_length = gene_coord[1] - gene_coord[0] + 1
                start_change = gene[0] - gene_coord[0] - 1
                incorrect_coords_dict[ordered_embl_dict[gene]] = {'length': prodigal_length,
                                                                  'start_change': start_change}
                break
    logger.debug('Number of genes with incorrect start predictions by Prodigal: ' +
                 str(len(list(incorrect_coords_dict.keys()))))
    pickle.dump(incorrect_coords_dict, open(dict_save_fp, "wb"))

    prodigal_correction_dict = incorrect_coords_dict
    ref_feature_dict = get_feature_by_locustag(ref_embl_record)
    global ref_sequence
    ref_sequence = ref_embl_record.seq
    ref_features = ref_embl_record.features

    global ref_prom_fp_dict
    ref_prom_fp_dict = get_prom_for_gene(ref_features, ref_sequence)

    for i in range(0, len(ratt_gbk_files)):
        ratt_contig_record = SeqIO.read(ratt_gbk_files[i], 'genbank')
        global record_sequence
        record_sequence = ratt_contig_record.seq
        prokka_noref_rec_pre = prokka_record_noref[i]
        prokka_noref_rec = correct_start_coords_prokka(prokka_noref_rec_pre, prodigal_correction_dict, record_sequence,
                                                       ref_sequence, ref_feature_dict, reference_locus_list,
                                                       reference_gene_locus_dict)
        global prokka_noref_dictionary
        prokka_noref_dictionary = generate_feature_dictionary(prokka_noref_rec.features)
        prokka_contig_record_pre = prokka_records[i]
        prokka_contig_record = correct_start_coords_prokka(prokka_contig_record_pre, prodigal_correction_dict,
                                                           record_sequence, ref_sequence, ref_feature_dict,
                                                           reference_locus_list, reference_gene_locus_dict)
        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        if i == 0 and not illumina:
            check_for_dnaA(ratt_contig_features)
        if illumina:
            logger.debug('DnaA might not be the first element. Circularization is not checked')
        global ratt_corrected_genes
        if len(ratt_correction_files) == 1:
            error_correction_fp = ratt_correction_files[0]
        else:
            try:
                error_correction_fp = ratt_correction_files[i]
                ratt_corrected_genes = get_ratt_corrected_genes(error_correction_fp, reference_gene_list)
            except IndexError:
                ratt_corrected_genes = []
        ratt_contig_non_cds = []
        for feature in ratt_contig_features:
            if feature.type == 'mRNA' or feature.type == 'rRNA':
                ratt_contig_non_cds.append(feature)
        logger.debug('Number of non-CDS elements')
        logger.debug(len(ratt_contig_non_cds))
        ratt_contig_features, ratt_blast_results = isolate_valid_ratt_annotations(ratt_contig_features,
                                                                                  ref_temp_fasta_dict,
                                                                                  reference_locus_list,
                                                                                  seq_ident=seq_ident,
                                                                                  seq_covg=seq_covg)
        ratt_contig_features = remove_duplicate_cds(ratt_contig_features)
        cds_from_ratt = 0
        for feature in ratt_contig_features:
            if feature.type == 'CDS':
                cds_from_ratt += 1
        logger.debug('Number of CDSs from RATT')
        logger.debug(str(cds_from_ratt))
        ratt_contig_features = get_ordered_features(ratt_contig_features)
        merged_genes, check_prokka = identify_merged_genes(ratt_contig_features)
        if check_prokka:
            merged_features, corner_cases, corner_cases_explicit, ratt_contig_features_mod, \
                prokka_contig_features_mod = get_annotation_for_merged_genes(merged_genes,
                                                                             prokka_contig_features,
                                                                             ratt_contig_features,
                                                                             reference_locus_list=reference_locus_list)
            ratt_contig_features = ratt_contig_features_mod
            prokka_contig_features = prokka_contig_features_mod
            merged_features_record = prokka_contig_record[:]
            merged_features_record.features = merged_features
            if corner_cases:
                logger.debug('MERGED GENES: Corner cases')
                for strand in corner_cases_explicit.keys():
                    if len(corner_cases_explicit[strand]) > 0:
                        logger.debug(isolate_id + ' '.join(corner_cases_explicit[strand]))
        else:
            merged_features = []
        if len(merged_features) > 0 and i == 0:
            SeqIO.write(merged_features_record, output_merged_genes, 'genbank')
        ratt_cds_count = 0
        for feature in ratt_contig_features:
            if feature.type == 'CDS':
                ratt_cds_count += 1
        global ratt_contig_features_dict
        ratt_contig_features_dict = generate_feature_dictionary(ratt_contig_features)
        global popped_ratt_genes
        popped_ratt_genes = []
        if len(ratt_contig_features) == 0:
            logger.warning("NO RATT ANNOTATION FOR CONTIG " + str(i + 1))
            feature_additions = {}
            feature_lengths = {}
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
            for prokka_feature in prokka_contig_record.features:
                if prokka_feature.type not in feature_additions.keys():
                    feature_additions[prokka_feature.type] = 1
                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                else:
                    feature_additions[prokka_feature.type] += 1
                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            continue
        elif len(prokka_contig_features) == 0:
            logger.warning("NO PROKKA ANNOTATION FOR CONTIG " + str(i + 1))
            prokka_contig_features = ratt_contig_features
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            logger.warning('Contig Number: ' + str(i + 1) + '\n')
            logger.warning('No Annotation to add from Prokka')
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the
            # features from RATT
            add_prokka_contig_record = prokka_contig_record[:]
            add_prokka_contig_record.features = []
            num_feat = 0
            add_prokka_contig_record.annotations['comment'] = 'Merged reference based annotation from RATT and ab ' \
                                                              'initio annotation from Prokka'
            ratt_annotation_mapping = {}  # Used for resolving annotations of overlapping features between RATT and
            # Prokka
            for index, feature in enumerate(ratt_contig_record.features):
                try:
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    ratt_annotation_mapping[(start, end)] = index
                except AttributeError:
                    logger.error('Attribute Error')
                    logger.error(feature)
                    logger.error(index)
            try:
                ratt_contig_record_mod = ratt_contig_record[:]
            except AttributeError:
                logger.error('Contains features with fuzzy locations')
                logger.error(ratt_contig_record)
            ratt_contig_record_mod.features = ratt_contig_features
            added_from_ratt = []
            non_cds_ratt = []
            global prokka_blast_list
            prokka_blast_list = []
            for feat in ratt_contig_record_mod.features:
                if feat.type == 'CDS':
                    added_from_ratt.append(feat.qualifiers['locus_tag'][0])
                else:
                    non_cds_ratt.append(feat)
            prokka_features_dict = generate_feature_dictionary(prokka_contig_features)
            prokka_features_not_in_ratt, ratt_overlapping_genes = \
                remove_duplicate_annotations(ratt_contig_features, prokka_features_dict)
            logger.debug('Number of prokka features not in RATT ' + str(len(list(prokka_features_not_in_ratt.keys()))))
            intergenic_ratt, intergenic_positions, ratt_pre_intergene, ratt_post_intergene = \
                get_interregions(ratt_contig_record_mod, intergene_length=1)
            sorted_intergenic_positions = sorted(intergenic_positions)
            feature_additions = {}
            feature_lengths = {}
            add_features_from_prokka = []
            for j in sorted_intergenic_positions:
                # Variable definitions
                ratt_unannotated_region_start = j[0]
                ratt_unannotated_region_end = j[1]
                ratt_strand = j[2]
                for feature_position in list(prokka_features_not_in_ratt.keys()):
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
                            mod_prokka_feature, invalid_prokka = \
                                validate_prokka_feature_annotation(prokka_feature, prokka_noref_dictionary,
                                                                   reference_gene_locus_dict, reference_locus_gene_dict,
                                                                   ref_temp_fasta_dict, ratt_blast_results,
                                                                   reference_locus_list=reference_locus_list,
                                                                   seq_ident=seq_ident,
                                                                   seq_covg=seq_covg)
                            add_features_from_prokka.append(mod_prokka_feature)
                        # If the Prokka feature overlaps with two RATT features
                        elif prokka_feature_start < ratt_unannotated_region_start and \
                                prokka_feature_end > ratt_unannotated_region_end:
                            if prokka_feature.type == 'source':  # This is to exclude the source feature as it is
                                # accounted for by RATT
                                break
                            ratt_overlapping_feature_1 = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                             prokka_strand)]
                            ratt_overlapping_feature_2 = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                              prokka_strand)]
                            ratt_overlapping_feature_1_loc = (int(ratt_overlapping_feature_1.location.start),
                                                              int(ratt_overlapping_feature_1.location.end),
                                                              int(ratt_overlapping_feature_1.location.strand))
                            ratt_overlapping_feature_2_loc = (int(ratt_overlapping_feature_2.location.start),
                                                              int(ratt_overlapping_feature_2.location.end),
                                                              int(ratt_overlapping_feature_2.location.strand))
                            locus_tag1 = ratt_overlapping_feature_1.qualifiers['locus_tag'][0]
                            locus_tag2 = ratt_overlapping_feature_2.qualifiers['locus_tag'][0]
                            overlapping_ratt_location_1 = {locus_tag1: ratt_overlapping_feature_1_loc}
                            overlapping_ratt_location_2 = {locus_tag2: ratt_overlapping_feature_2_loc}
                            overlapping_ratt_locus_tags = [locus_tag1, locus_tag2]
                            mod_prokka_feature, invalid_prokka = \
                                validate_prokka_feature_annotation(prokka_feature, prokka_noref_dictionary,
                                                                   reference_gene_locus_dict, reference_locus_gene_dict,
                                                                   ref_temp_fasta_dict, ratt_blast_results,
                                                                   seq_ident=seq_ident,
                                                                   seq_covg=seq_covg,
                                                                   check_ratt=True,
                                                                   ratt_features=overlapping_ratt_locus_tags[0],
                                                                   ratt_locations=overlapping_ratt_location_1,
                                                                   reference_locus_list=reference_locus_list)
                            if not invalid_prokka:
                                mod_prokka_feature, invalid_prokka_2 = \
                                    validate_prokka_feature_annotation(prokka_feature, prokka_noref_dictionary,
                                                                       reference_gene_locus_dict,
                                                                       reference_locus_gene_dict,
                                                                       ref_temp_fasta_dict,
                                                                       ratt_blast_results,
                                                                       seq_ident=seq_ident,
                                                                       seq_covg=seq_covg,
                                                                       check_ratt=True,
                                                                       ratt_features=overlapping_ratt_locus_tags[1],
                                                                       ratt_locations=overlapping_ratt_location_2,
                                                                       reference_locus_list=reference_locus_list)
                                if not invalid_prokka_2:
                                    add_features_from_prokka.append(mod_prokka_feature)
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
                            does_not_overlap = False
                            if (prokka_feature_start < ratt_unannotated_region_start) and \
                                    (prokka_feature_end in ratt_unannotated_region_range):
                                ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                               prokka_strand)]
                            elif (prokka_feature_start in ratt_unannotated_region_range) and \
                                    prokka_feature_end > ratt_unannotated_region_end:
                                ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                                prokka_strand)]
                            else:
                                does_not_overlap = True
                            if not does_not_overlap:
                                overlapping_ratt_locus_tags = [ratt_overlapping_feature.qualifiers['locus_tag'][0]]
                                ratt_overlapping_feature_loc = (int(ratt_overlapping_feature.location.start),
                                                                int(ratt_overlapping_feature.location.end),
                                                                int(ratt_overlapping_feature.location.strand))
                                overlapping_ratt_locations = {ratt_overlapping_feature.qualifiers['locus_tag'][0]:
                                                              ratt_overlapping_feature_loc}
                                # If Prokka feature is annotated as an L_tag and not a gene, blast it with
                                # overlapping Rv and then Blast it with reference to get the closer hit
                                if prokka_feature.type == 'CDS' and 'gene' not in prokka_feature.qualifiers.keys():
                                    check_if_same_gene = blast_ratt_anno_with_prokka_ltag(ratt_overlapping_feature,
                                                                                          prokka_feature,
                                                                                          seq_ident=seq_ident,
                                                                                          seq_covg=seq_covg)
                                    if check_if_same_gene:
                                        prokka_feature.qualifiers['gene'] = \
                                            ratt_overlapping_feature.qualifiers['locus_tag']

                                add_prokka_contig_record, \
                                    included = check_inclusion_criteria(ratt_annotation_mapping,
                                                                        add_prokka_contig_record,
                                                                        ratt_overlapping_feature,
                                                                        prokka_feature,
                                                                        overlapping_ratt_locus_tags,
                                                                        overlapping_ratt_locations,
                                                                        ref_temp_fasta_dict=ref_temp_fasta_dict,
                                                                        reference_gene_locus_dict=reference_gene_locus_dict,
                                                                        reference_locus_gene_dict=reference_locus_gene_dict,
                                                                        ratt_blast_results=ratt_blast_results,
                                                                        reference_locus_list=reference_locus_list,
                                                                        seq_ident=seq_ident,
                                                                        seq_covg=seq_covg)
                                if included:  # To check if Prokka feature passed the inclusion criteria and was
                                    # integrated into the EMBL file.
                                    # The if-else condition below is to keep track of the features added from Prokka
                                    # for the log file
                                    if prokka_feature.type not in feature_additions.keys():
                                        feature_additions[prokka_feature.type] = 1
                                        feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                                    else:
                                        feature_additions[prokka_feature.type] += 1
                                        feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            if len(merged_features) > 0:
                for feature_1 in merged_features:
                    add_prokka_contig_record.features.append(feature_1)
            for prokka_feature_append in add_features_from_prokka:
                add_prokka_contig_record.features.append(prokka_feature_append)
            for ratt_feature_append in ratt_contig_features_dict.values():
                add_prokka_contig_record.features.append(ratt_feature_append)
            for non_cds in ratt_contig_non_cds:
                add_prokka_contig_record.features.append(non_cds)
            records_ref, positions_ref, pre_ref, post_ref = get_interregions(add_prokka_contig_record,
                                                                             intergene_length=1)
            add_features_from_prokka_ref = []
            sorted_positions_ref = sorted(positions_ref)
            prokka_ref_features_sorted = get_ordered_features(prokka_contig_features)
            for fill_pos in sorted_positions_ref:
                if fill_pos[2] == '+':
                    fill_pos_strand = 1
                else:
                    fill_pos_strand = -1
                for ref_feat in prokka_ref_features_sorted:
                    ref_feat_start = int(ref_feat.location.start)
                    ref_feat_end = int(ref_feat.location.end)
                    ref_feat_strand = int(ref_feat.location.strand)
                    if ref_feat_strand != fill_pos_strand:
                        continue
                    elif ref_feat.type != 'CDS':
                        continue
                    elif ref_feat_start > fill_pos[1] and ref_feat_end > fill_pos[1]:
                        break
                    elif ref_feat_start < fill_pos[0] and ref_feat_end < fill_pos[0]:
                        continue
                    elif ref_feat_start > fill_pos[0] and ref_feat_end < fill_pos[1]:
                        if ref_feat.qualifiers['product'][0] == 'hypothetical protein':
                            continue
                        if 'gene' in ref_feat.qualifiers.keys() and '_' in ref_feat.qualifiers['gene'][0]:
                            ref_feat.qualifiers['gene'] = ref_feat.qualifiers['locus_tag']
                        elif 'gene' not in ref_feat.qualifiers.keys():
                            ref_feat.qualifiers['gene'] = ref_feat.qualifiers['locus_tag']
                        else:
                            mod_ref_feat, \
                            invalid_ref_feat = validate_prokka_feature_annotation(ref_feat,
                                                                                  prokka_noref_dictionary,
                                                                                  reference_gene_locus_dict,
                                                                                  reference_locus_gene_dict,
                                                                                  ref_temp_fasta_dict,
                                                                                  ratt_blast_results,
                                                                                  reference_locus_list,
                                                                                  seq_ident=seq_ident,
                                                                                  seq_covg=seq_covg)
                            ref_feat = mod_ref_feat
                        add_ref_note = 'This annotation is added from Prokka reference run'
                        if 'note' in ref_feat.qualifiers.keys():
                            ref_feat.qualifiers['note'].append(add_ref_note)
                        else:
                            ref_feat.qualifiers['note'] = [add_ref_note]
                        add_features_from_prokka_ref.append(ref_feat)

            for prokka_ref_feat in add_features_from_prokka_ref:
                add_prokka_contig_record.features.append(prokka_ref_feat)
            annomerge_records.append(add_prokka_contig_record)
    # Post-processing of genbank file to remove duplicates and rename locus_tag for
    # Prokka annotations
    prokka_record_fp = file_path + 'prokka-noreference/' + isolate_id + '.gbk'
    prokka_record_noref = list(SeqIO.parse(prokka_record_fp, 'genbank'))
    annomerge_records_post_processed = []
    for rec_num in range(0, len(annomerge_records)):
        prokka_rec = annomerge_records[rec_num]
        prokka_noref_rec_pre = prokka_record_noref[rec_num]
        prokka_noref_rec = correct_start_coords_prokka(prokka_noref_rec_pre, prodigal_correction_dict, record_sequence,
                                                       ref_sequence, ref_feature_dict, reference_locus_list,
                                                       reference_gene_locus_dict)
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
        checked_features = 0
        for feature in raw_features:
            checked_features += 1
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
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    if loc_key in added_cds.keys():
                        continue
                    if loc_key in prokka_noref_dict.keys():
                        feature_noref = prokka_noref_dict[loc_key]
                        if 'gene' in feature_noref.qualifiers.keys():
                            if 'gene_synonym' in feature_noref.qualifiers.keys():
                                feature_noref.qualifiers['gene_synonym'].append(feature_noref.qualifiers['gene'])
                            else:
                                feature_noref.qualifiers['gene_synonym'] = feature_noref.qualifiers['gene']
                            feature_noref.qualifiers['gene'] = feature_noref.qualifiers['locus_tag']
                    else:
                        feature_noref = feature
                    if 'gene' in feature_noref.qualifiers.keys() and '_' in feature_noref.qualifiers['gene'][0]:
                        new_locus_tag = 'L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['locus_tag'] = [new_locus_tag]
                        new_protein_id = 'C:L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['protein_id'] = [new_protein_id]
                        feature_noref.qualifiers['gene'] = feature_noref.qualifiers['locus_tag']
                    else:
                        new_locus_tag = 'L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['locus_tag'] = [new_locus_tag]
                        new_protein_id = 'C:L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                        feature_noref.qualifiers['protein_id'] = [new_protein_id]
                    added_cds[loc_key] = feature_noref
                    prokka_rec.features.append(feature_noref)
                else:
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    if loc_key in added_cds.keys():
                        continue
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
        # Blasting unannotated CDSs from Prokka with H37Rv
        # TODO: Delete this part. This is for internal use where unannotated genes are blasted to previously assigned
        # novel genes to transfer annotation.
        if check_mtb and len(mtb_fasta_fp) == 0:
            logger.warning('File path for novel genes not specified')
        elif check_mtb:
            mtb_fasta = mtb_fasta_fp
            for feature in prokka_rec.features:
                if feature.type == 'CDS' and (('gene' not in feature.qualifiers.keys() and
                                               ('L_' in feature.qualifiers['locus_tag'][0] or
                                                'L2_' in feature.qualifiers['locus_tag'][0])) or
                                              ('gene' in feature.qualifiers.keys() and
                                               '_' in feature.qualifiers['gene'][0] and
                                               'L2_' in feature.qualifiers['locus_tag'][0])):
                    query_sequence = feature.qualifiers['translation'][0]
                    blast_to_mtb = NcbiblastpCommandline(subject=mtb_fasta, outfmt='"7 qseqid qlen sseqid slen qlen '
                                                                                   'length pident qcovs"')
                    stdout_2, stderr_2 = blast_to_mtb(stdin=query_sequence)
                    mtb_hit, all_hits_mtb = identify_top_hits_mtb(stdout_2, identity=seq_ident, coverage=seq_covg, mtb=True)
                    for gene in all_hits_mtb.keys():
                        note = 'locus_tag:' + gene + ':' + str(all_hits_mtb[gene])
                        if 'note' in feature.qualifiers.keys():
                            feature.qualifiers['note'].append(note)
                        else:
                            feature.qualifiers['note'] = [note]
        # Adding translated sequence to RATT annotations
        for feature in prokka_rec.features:
            if feature.type == 'CDS' and 'translation' not in feature.qualifiers.keys():
                feature_sequence = translate(feature.extract(record_sequence), table=11, to_stop=True)
                feature.qualifiers['translation'] = [feature_sequence]
        # Verifying annotations for essential genes
        if essentiality:
            essential_genes = get_essential_genes()
            all_rv_genes_in_isolate = []
            for feature in prokka_rec.features:
                if feature.type == 'CDS':
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    if locus_tag[:2] == 'Rv':
                        all_rv_genes_in_isolate.append(locus_tag)
                        if 'gene' not in feature.qualifiers.keys() and locus_tag in reference_locus_gene_dict.keys():
                            feature.qualifiers['gene'] = [reference_locus_gene_dict[locus_tag]]
                    new_locus_tag = feature.qualifiers['locus_tag'][0]
                    if 'gene' not in feature.qualifiers.keys():
                        feature.qualifiers['gene'] = [new_locus_tag]
                    else:
                        if '_' in feature.qualifiers['gene'][0] and (feature.qualifiers['gene'][0][:2] != 'L_' or
                                                                 feature.qualifiers['gene'][0][:3] != 'L2_'):
                            if '_' not in feature.qualifiers['locus_tag'][0]:
                                feature.qualifiers['locus_tag'] = [new_locus_tag]
                            else:
                                new_locus_tag = 'L2_' + feature.qualifiers['locus_tag'][0].split('_')[1]
                                feature.qualifiers['locus_tag'] = [new_locus_tag]
                                feature.qualifiers['gene'] = [new_locus_tag]
                    if new_locus_tag in essential_genes:
                        rv_length = ref_protein_lengths[new_locus_tag]
                        feature_length = len(feature.qualifiers['translation'][0])
                        feature_coverage = float(feature_length)/float(rv_length)
                        if feature_coverage <= 0.95:
                            logger.warning('Essential gene ' + new_locus_tag + ' is truncated in isolate ' +
                                         isolate_id + '. Coverage: ' + str(feature_coverage*100))
        # Get remaining unannotated region
        if add_noref_annotations:
            ordered_features_final = get_ordered_features(prokka_rec.features)
            prokka_rec.features = ordered_features_final
            records, positions, pre, post = get_interregions(prokka_rec, intergene_length=1)
            positions_lengths = [len(p) for p in records]
            add_features_from_prokka_noref = []
            sorted_positions = sorted(positions)
            prokka_noref_features_sorted = get_ordered_features(prokka_noref_rec.features)
            for fill_pos in sorted_positions:
                if fill_pos[2] == '+':
                    fill_pos_strand = 1
                else:
                    fill_pos_strand = -1
                for noref_feat in prokka_noref_features_sorted:
                    noref_feat_start = int(noref_feat.location.start)
                    noref_feat_end = int(noref_feat.location.end)
                    noref_feat_strand = int(noref_feat.location.strand)
                    if noref_feat_strand != fill_pos_strand:
                        continue
                    elif noref_feat.type != 'CDS':
                        continue
                    elif noref_feat_start > fill_pos[1] and noref_feat_end > fill_pos[1]:
                        break
                    elif noref_feat_start < fill_pos[0] and noref_feat_end < fill_pos[0]:
                        continue
                    elif noref_feat_start > fill_pos[0] and noref_feat_end < fill_pos[1]:
                        if noref_feat.qualifiers['product'][0] == 'hypothetical protein':
                            continue
                        if 'gene' in noref_feat.qualifiers.keys() and '_' in noref_feat.qualifiers['gene'][0]:
                            noref_feat.qualifiers['gene'] = noref_feat.qualifiers['locus_tag']
                        elif 'gene' not in noref_feat.qualifiers.keys():
                            noref_feat.qualifiers['gene'] = noref_feat.qualifiers['locus_tag']
                        add_noref_note = 'This annotation is added from Prokka no-reference run'
                        if 'note' in noref_feat.qualifiers.keys():
                            noref_feat.qualifiers['note'].append(add_noref_note)
                        else:
                            noref_feat.qualifiers['note'] = [add_noref_note]
                        add_features_from_prokka_noref.append(noref_feat)
                        prokka_rec.features.append(noref_feat)
            logger.debug('To add from noref: ' + str(len(add_features_from_prokka_noref)))

        # Removing gene names assigned based on domains
        logger.debug('Final feature annotation verification')
        added_ltags = []
        for feature_final in prokka_rec.features:
            if feature_final.type != 'CDS':
                continue
            if 'inference' in feature_final.qualifiers.keys():
                is_domain = False
                for infer in feature_final.qualifiers['inference']:
                    if infer.startswith('protein motif'):
                        is_domain = True
                        break
                if is_domain:
                    feature_final.qualifiers['gene'] = feature_final.qualifiers['locus_tag']
            # Check if same locus tag is annotated with 2 different sequences
            final_feature_name = feature_final.qualifiers['locus_tag'][0]
            if final_feature_name in added_ltags and (final_feature_name.startswith('L_') or
                                                      final_feature_name.startswith('L2_')):
                prev_ltag = final_feature_name.split('_')[1]
                ltag = final_feature_name.split('_')[0][-1]
                if ltag == 'L':
                    renamed_ltag = 'L2_' + str(prev_ltag)
                else:
                    add_digit = int(ltag) + 1
                    renamed_ltag = 'L' + str(add_digit) + '_' + str(prev_ltag)
                feature_final.qualifiers['locus_tag'] = [renamed_ltag]
            if 'gene' not in feature_final.qualifiers.keys():
                feature_final.qualifiers['gene'] = feature_final.qualifiers['locus_tag']
            added_ltags.append(feature_final.qualifiers['locus_tag'][0])
        sorted_final = get_ordered_features(prokka_rec.features)
        prokka_rec.features = sorted_final
        annomerge_records_post_processed.append(prokka_rec)
        output_isolate_recs = [r for r in annomerge_records_post_processed]
        isolate_features = annomerge_records_post_processed[0].features[:]
        output_isolate_recs[0].features = []
        prev_feature_list = []
        num_overlaps = 0
        positions_to_be_resolved = []
        resolve_pairs = []
        for feature in isolate_features:
            if feature.type != 'CDS':
                continue
            if len(prev_feature_list) == 0:
                prev_feature_list = [int(feature.location.start), int(feature.location.end),
                                     int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                prev_feature = feature
                continue
            if int(feature.location.start) <= prev_feature_list[1] and \
                    int(feature.location.start) >= prev_feature_list[0] and \
                    (int(feature.location.start) == prev_feature_list[0] or
                     int(feature.location.end) == prev_feature_list[1]):
                num_overlaps += 1
                positions_to_be_resolved.append((prev_feature_list[0], prev_feature_list[1], prev_feature_list[2]))
                prev_feature_list = [int(feature.location.start), int(feature.location.end),
                                     int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                positions_to_be_resolved.append(
                    (int(feature.location.start), int(feature.location.end), int(feature.location.strand)))
                resolve_pairs.append([prev_feature, feature])
                prev_feature = feature
            elif int(feature.location.end) <= prev_feature_list[1] and \
                    int(feature.location.end) >= prev_feature_list[0] and \
                    (int(feature.location.start) == prev_feature_list[0] or
                     int(feature.location.end) == prev_feature_list[1]):
                num_overlaps += 1
                positions_to_be_resolved.append((prev_feature_list[0], prev_feature_list[1], prev_feature_list[2]))
                prev_feature_list = [int(feature.location.start), int(feature.location.end),
                                     int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                positions_to_be_resolved.append(
                    (int(feature.location.start), int(feature.location.end), int(feature.location.strand)))
                resolve_pairs.append([prev_feature, feature])
                prev_feature = feature
            else:
                prev_feature_list = [int(feature.location.start), int(feature.location.end),
                                     int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                prev_feature = feature
        for feature in isolate_features:
            if feature.type != 'CDS':
                output_isolate_recs[0].features.append(feature)
            else:
                position = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                if position not in positions_to_be_resolved:
                    output_isolate_recs[0].features.append(feature)
        for feat_pair in resolve_pairs:
            if 'protein_id' in feat_pair[0].qualifiers.keys():
                prokka_annotation = feat_pair[0]
                ratt_annotation = feat_pair[1]
            else:
                prokka_annotation = feat_pair[1]
                ratt_annotation = feat_pair[0]
            take_ratt = pick_best_hit(ratt_annotation, prokka_annotation, isolate_sequence)
            if take_ratt:
                output_isolate_recs[0].features.append(ratt_annotation)
            else:
                output_isolate_recs[0].features.append(prokka_annotation)
        ordered_feats = get_ordered_features(output_isolate_recs[0].features)
        output_isolate_recs[0].features = ordered_feats[:]

    SeqIO.write(output_isolate_recs, output_genbank, 'genbank')
    final_cdss = [f for f in output_isolate_recs[0].features if f.type == 'CDS']
    logger.debug('Number of CDSs annomerge: ' + str(len(final_cdss)))
    logger.debug('annomerge run time: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes')
