__author__ = "Deepika Gunasekaran"
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
import functools
import collections
import os
import tempfile
import pickle
import logging
import time
import re
import subprocess
from copy import deepcopy
from math import log, ceil

# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition, SeqFeature
from Bio import Align

from . import BLAST
from . import converter
from . import config
from . import designator
from . import extractor
from . import __version__


def overlap_inframe(loc1, loc2):
    """
    Say whether two FeatureLocations are overlapping and share the same reading frame.
    This is for the purpose of CDSs to determine whether they should be corresponding to the
    same genomic feature or whether one gene has lost its stop codon and merged into its neighbor.

    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features overlap in-frame
    """
    # overlap determination adapted from https://stackoverflow.com/a/2953979
    overlap = (min(loc1.end, loc2.end) - max(loc1.start, loc2.start))

    if loc1.strand == loc2.strand and overlap > 0:
        # pseudogenes may occupy multiple reading frames, so
        # check both the start-defined and stop-defined frames,
        # as well as the actual overlapping part.
        if (((loc1.start - loc2.start) % 3 == 0
             or (loc1.end - loc2.end) % 3 == 0)
            and (overlap % 3 == 0)
        ):
            return True
    return False

def have_same_stop(loc1, loc2):
    """
    Say whether two FeatureLocations have the same stop position.
    This is a special case of overlap_inframe().
    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features have the same stop position
    """
    if loc1.strand == loc2.strand:
        return (
            (loc1.strand == 1 and loc1.end == loc2.end)
            or
            (loc1.strand == -1 and loc1.start == loc2.start)
        )
    return False

def is_broken_stop(feature):
    """
    Finds the amount and location of internal stops.
    :param feature: A SeqFeature object
    """
    internal_stop = False
    note = ''
    translation = str(feature.extract(record_sequence).translate(to_stop=False))
    num_stop = [i for i,e in enumerate(translation) if e == "*"]
    if len(num_stop) > 1 or translation[-1] != "*":
        internal_stop = True
        note = f"Hybran: Internal stop detected in the following codon(s): {','.join([str(i) for i in num_stop])}"
    return internal_stop, note

def log_feature_fate(feature, logfile, remark=""):
    """
    General-purpose logging function to print out a gene's information and a comment
    :param feature: A SeqFeature object
    :param logfile: An open filehandle
    :param remark: (str) A comment
    """
    if 'locus_tag' in feature.qualifiers:
        locus_tag = feature.qualifiers['locus_tag'][0]
    else:
        locus_tag = feature.id
    print('\t'.join([locus_tag, remark]), file=logfile)

def log_coord_correction(og_feature, feature, logfile):
    """
    This function is used to log gene information when coord_check() fixes a start/stop postition
    :param og_feature: Original SeqFeature object
    :param feature: Updated SeqFeature object
    :param logfile: An open filehandle
    """
    locus_tag = og_feature.qualifiers['locus_tag'][0]
    gene_name = og_feature.qualifiers['gene'][0]
    strand = str(og_feature.strand)
    og_start = str(int(og_feature.location.start) + 1)
    og_end = (og_feature.location.end)
    new_start = str(int(feature.location.start) + 1)
    new_end = (feature.location.end)
    start_fixed = str(og_start != new_start).lower()
    stop_fixed = str(og_end != new_end).lower()

    if designator.is_pseudo(og_feature.qualifiers) == designator.is_pseudo(feature.qualifiers):
        if designator.is_pseudo(og_feature.qualifiers):
            status = "remains_pseudo"
        else:
            status = "remains_non_pseudo"
    elif designator.is_pseudo(feature.qualifiers):
        status = "became_pseudo"
    else:
        status = "became_non_pseudo"

    if og_feature.strand == -1:
        start_fixed, stop_fixed = stop_fixed, start_fixed
    line = [locus_tag, gene_name, strand, og_start, og_end, new_start, new_end,  start_fixed, stop_fixed, status]
    print('\t'.join(str(v) for v in line), file=logfile)


def load_reference_info(proteome_fasta):
    """
    This function gets all the global variables for annomerge.
    1. A list of genes from the reference protein fasta
    2. A list of locus tags from the reference protein fasta
    3. A dict with gene name as keys and locus tags as values
    4. A dict with locus tags as keys and gene name as values
    5. A dict with locus tags as keys and amino acid sequence as values
    6. A dict with locus tags as keys and amino acid lengths as values
    7. A dict with locus tags as keys and nucleotide lengths as values
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    reference_gene_list = []
    reference_locus_list = []
    reference_locus_gene_dict = {}
    reference_gene_locus_dict = collections.defaultdict(list)
    ref_protein_lengths = {}

    ref_fasta_records = SeqIO.parse(proteome_fasta, 'fasta')
    for record in ref_fasta_records:
        gene_locus_name = record.id
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
        reference_gene_locus_dict[gene].append(locus)
        reference_locus_gene_dict[locus] = gene
        seq = str(record.seq)
        ref_protein_lengths[locus] = len(seq)
    return reference_gene_list, reference_locus_list, reference_gene_locus_dict, reference_locus_gene_dict, ref_protein_lengths


def upstream_context(feature_location, source_seq, n=40, circular=True):
    """
    Get the n bases upstream of a genomic feature.
    :param feature_location: FeatureLocation
    :source_seq: Seq or str sequence to which the FeatureLocation maps
    :param n: int number of bases of upstream context to retrieve
    :param circular: boolean whether the end of the sequence precedes position 0
                     (linear support not currently implemented)
    :return: str upstream sequence length n
    """
    feature_start = int(feature_location.start)
    feature_stop = int(feature_location.end)
    feature_strand = int(feature_location.strand)
    if feature_strand == 1:
        prom_end = feature_start
        if feature_start < n:
            if circular:
                prev_prom_len = n - prom_end
                prom_start = len(source_seq) - prev_prom_len
                extra_prom_seq = str(source_seq[prom_start:len(source_seq)])
            else:
                extra_prom_seq = ''
            prom_seq = extra_prom_seq + str(source_seq[0:prom_end])
        else:
            prom_start = feature_start - n
            prom_seq = str(source_seq[prom_start:prom_end])
    else:
        prom_start = feature_stop
        if len(source_seq) - feature_stop < n:
            if circular:
                prev_prom_len = len(source_seq) - prom_start
                extra_prom_seq = str(source_seq[0:n-prev_prom_len])
            else:
                extra_prom_seq = ''
            prom_seq = str(source_seq[prom_start:]) + extra_prom_seq
        else:
            prom_end = feature_stop + n
            prom_seq = str(source_seq[prom_start:prom_end])
        prom_seq = str(Seq(prom_seq).reverse_complement())

    return prom_seq

def get_nuc_seq_for_gene(feature_list, source_seq):
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
    ref_fna_dict = {}
    for f in feature_list:
        if f.type != 'CDS':
            continue
        locus = f.qualifiers['locus_tag'][0]
        prom_seq = upstream_context(f.location, source_seq)
        ref_fna_dict[locus] = SeqRecord(seq=f.extract(source_seq), id=locus)
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
    return ref_prom_dict, ref_fna_dict


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
        # If the assigned new gene name is in the reference, increment the alphabet suffix
        char_start += 1
        if strand == '-1' or strand == '-' or strand == -1:
            new_gene_name = gene_name + chr(char_start) + 'c'
        else:
            new_gene_name = gene_name + chr(char_start)
    return new_gene_name

def merge_qualifiers(f1quals, f2quals):
    """
    Combine two qualifier dictionaries, combining lists for qualifiers that may
    have multiple values.
    :param f1quals: dict first SeqFeature's qualifiers
    :param f2quals: dict second SeqFeature's qualifiers
    :return: dict combined qualifier dictionary
    """
    multifields = [
            'note',
            'gene_synonym',
            'experiment',
            'inference',
    ]
    final_qualifiers = deepcopy(f1quals)
    final_qualifiers.update(f2quals)
    for qual in multifields:
        if qual in f1quals.keys() and qual in f2quals.keys():
            final_qualifiers[qual] = list(
                set(f1quals[qual]).union(set(f2quals[qual]))
            )
    return final_qualifiers

def process_split_genes(flist):
    """
    Given a list of features ordered by genomic position, assign the same
    locus tag to consecutive fragments of the same gene.
    :param flist: list of SeqFeature objects
    :returns:
        list of SeqFeature objects to keep (some modified from the original)
        list of annotations that have been merged into their neighbor.
    """
    logger = logging.getLogger('ProcessSplitGenes')
    outlist = []
    dropped_ltag_features = []
    last_gene_by_strand = {}
    for feature in flist:
        outlist.append(feature)
        if feature.type != 'CDS':
            continue
        if feature.location.strand in last_gene_by_strand:
            last_gene = last_gene_by_strand[feature.location.strand]
            # we'll be working with the loop variable so we can update the dictionary ahead of time
            last_gene_by_strand[feature.location.strand] = feature
        else:
            last_gene_by_strand[feature.location.strand] = feature
            continue

        last_gene_named = 'gene' in last_gene.qualifiers
        curr_gene_named = 'gene' in feature.qualifiers
        only_one_named = last_gene_named ^ curr_gene_named

        combine = False
        reason = ''

        # Many ab initio CDSs are rejected before this step, so we cannot deduce adjacency on the genome by adjacency on the list.
        # We can use the difference in their locus tag numbers to determine that.
        n_last_gene = int(last_gene.qualifiers['locus_tag'][0].split('_')[1])
        n_curr_gene = int(feature.qualifiers['locus_tag'][0].split('_')[1])
        dist = n_curr_gene - n_last_gene

        #
        # This function runs after coordinate correction, so if two ab initio CDSs now overlap in-frame,
        # it means that one of them had its start position corrected to correspond to that of the other fragment.
        #
        if overlap_inframe(last_gene.location, feature.location):
            combine = True
            reason = 'overlapping_inframe'
        #
        # Don't merge CDSs that are too far apart.
        # We only want to consider those that are directly adjacent, but we have to account for alternating genes on the opposite strand.
        # Using 3 as a cutoff because I do not expect more than one or two overlapping genes/gene fragments on the opposite strand between two fragments of what should be a single gene.
        #
        elif dist > 3:
            continue
        #
        # When neither are named, we have no other way to know whether they should be combined
        #
        elif (not last_gene_named) and (not curr_gene_named):
            continue
        #
        # Check for complementarity
        #
        elif only_one_named or extractor.get_gene(last_gene) == extractor.get_gene(feature):
            if only_one_named and last_gene_named:
                ref_gene = last_gene.qualifiers['gene'][0]
            else:
                ref_gene = feature.qualifiers['gene'][0]

            lg_status = coord_check(
                last_gene,
                ref_annotation[ref_gene],
            )
            cg_status = coord_check(
                feature,
                ref_annotation[ref_gene],
            )
            if ((any(lg_status) and any(cg_status))
                and (int(lg_status[0])+int(cg_status[0]), int(lg_status[1])+int(cg_status[1]))==(1,1)):
                combine = True
                reason = 'complementary_fragments'

        #
        # Combine intervals and check validity, aborting otherwise
        #
        if combine:
            new_start = last_gene.location.start
            new_end = feature.location.end
            if only_one_named and last_gene_named:
                new_feature = deepcopy(last_gene)
                dropped_feature = feature
            else:
                new_feature = deepcopy(feature)
                dropped_feature = last_gene
            dropped_feature_name = f"{extractor.get_ltag(dropped_feature)}:{extractor.get_gene(dropped_feature)}"
            new_feature_name = f"{extractor.get_ltag(new_feature)}:{extractor.get_gene(new_feature)}"
            new_feature.location = FeatureLocation(new_start, new_end, feature.location.strand)
            new_feature.qualifiers = merge_qualifiers(dropped_feature.qualifiers, new_feature.qualifiers)
            new_feature.qualifiers['pseudo'] = ['']

            if all(coord_check(new_feature, ref_annotation[new_feature.qualifiers['gene'][0]], fix_start=True, fix_stop=True)) or reason == 'overlapping_inframe':
                dropped_ltag_features.append(
                    (dropped_feature, f"{dropped_feature_name} combined with {new_feature_name}: {reason}")
                )
                last_gene_by_strand[feature.location.strand] = new_feature
                outlist.remove(last_gene)
                outlist.remove(feature)
                outlist.append(new_feature)
            else:
                logger.debug(f"Did not combine {dropped_feature_name} and {new_feature_name} ({reason}) due to failing coordinate verification.")

    return outlist, dropped_ltag_features

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
    merged_gene_locus = collections.defaultdict(list)
    merged_features_addition = []
    features_from_ratt = collections.defaultdict(list)
    genes_from_ratt = collections.defaultdict(list)
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
        elif ratt_feature_location in check_positions_for_annotation[feature.strand]:
            merged_gene_locus[ratt_feature_location].append('|'.join([extractor.get_ltag(feature), extractor.get_gene(feature)]))
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
                designator.append_qualifier(
                    feature.qualifiers, 'note',
                    'Fusion: ' + merged_features_string
                )
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
                features_from_ratt[feature_location].append(feature)
                genes_from_ratt[feature_location].append('|'.join([extractor.get_ltag(feature), extractor.get_gene(feature)]))
    for location in features_from_ratt.keys():
        new_feature = features_from_ratt[location][0]
        merged_features_string = ",".join(genes_from_ratt[location])
        new_gene_name = rename_locus(new_feature.qualifiers['locus_tag'][0],
                                     new_feature.location.strand,
                                     reference_locus_list)
        new_feature.qualifiers['locus_tag'] = [new_gene_name]
        designator.append_qualifier(
            new_feature.qualifiers, 'note',
            'Fusion: ' + merged_features_string
        )
        merged_features_addition.append(new_feature)
        merged_genes[new_feature.location.strand].remove(location)
    if len(merged_genes[-1]) > 0 or len(merged_genes[1]) > 0:
        corner_cases = True
    else:
        corner_cases = False
    return merged_features_addition, corner_cases, merged_genes, final_ratt_features, final_prokka_features


def identify_conjoined_genes(ratt_features):
    """
    This function identifies RATT annotations for genes that had a nonstop mutation and combined with the following gene, adds a note and tags the feature 'pseudo'.
    It differs from identify_merged_genes() in that the latter looks for RATT artifacts where the result is two annotations with the same coordinates, while this function does not mark any annotations for removal.

    :param ratt_features: list of sorted SeqFeature objects. This object will be modified.
    :return: list of the SeqFeature objects that were identified and tagged.
    """
    prev_feature = None
    conjoined_features = []
    for feature in ratt_features:
        if not prev_feature:
            prev_feature = feature
            continue

        if (prev_feature.location != feature.location # dealing with these is the job of identify_merged_genes()
            and have_same_stop(prev_feature.location, feature.location)
        ):
            if feature.strand == -1:
                upstream = feature
                downstream = prev_feature
            else:
                upstream = prev_feature
                downstream = feature

            designator.append_qualifier(
                upstream.qualifiers,
                'note',
                f"Nonstop mutation in this gene causes an in-frame fusion with {'|'.join([extractor.get_ltag(downstream),extractor.get_gene(downstream)])}"
            )
            upstream.qualifiers['pseudo'] = ['']
            conjoined_features.append(upstream)

        prev_feature = feature

    return conjoined_features


def liftover_annotation(feature, ref_feature, pseudo, inference):
    """
    Add ref_feature's functional annotation to feature.

    :param feature: SeqFeature ab initio annotation.
                    This argument is modified by this function.
    :param ref_feature: SeqFeature reference annotation
    :param pseudo: bool whether feature should have the `pseudo` qualifier
    :param inference: str /inference annotation justifying the liftover
    """

    for rubbish in ['gene', 'protein_id']:
        feature.qualifiers.pop(rubbish, None)
    # Remove inferences for the assignments that we're going to discard.
    # Only keep the ab initio inference from the ORF finder and ones that we assign.
    feature.qualifiers['inference'][:] = [
        _ for _ in feature.qualifiers['inference']
        if ('ab initio prediction' in _ or
            'Hybran' in _
            )
    ]
    # Add our own qualifier
    feature.qualifiers['inference'].append(
        inference
    )

    # make a copy of the reference qualifiers dict so
    # we can modify it before merging (we don't want to carry
    # over certain attributes)
    ref_feature_qualifiers_copy = deepcopy(ref_feature.qualifiers)
    ref_specific = [
        'locus_tag',
        'old_locus_tag',
        'translation',
    ]
    for qual in ref_specific:
        ref_feature_qualifiers_copy.pop(qual, None)

    feature.qualifiers = merge_qualifiers(
        feature.qualifiers,
        ref_feature_qualifiers_copy,
    )

    # avoid inheriting a pseudo tag from the reference if
    # it's not warranted
    if pseudo:
        feature.qualifiers['pseudo'] = ['']
        feature.qualifiers.pop('translation', None)
        designator.append_qualifier(
            feature.qualifiers, 'note',
            'Pseudo: Low subject coverage, but passing query coverage.')
    else:
        feature.qualifiers.pop('pseudo', None)
        feature.qualifiers.pop('pseudogene', None)

def coord_check(feature, ref_feature, fix_start=False, fix_stop=False, ref_gene_name=None
):
    """
    This function takes a feature as an input and aligns it to the corresponding reference gene.
    This function will return two Boolean values (True if the start/end coordinates align with the reference).
    The start/end coordinates can be corrected to match the reference if desired in various circumstances.
    :param feature: SeqFeature object
    :param fix_start: Boolean
    :param fix_stop: Boolean
    :param ref_gene_name: str reference gene to check against.
        Must be a key existing in the `ref_annotation` dictionary.
        If not defined, the reference gene matching `feature`'s gene qualifier is used instead.
    :return: True/False if the start/stop was fixed
    """

    ref_seq = ref_feature.extract(ref_sequence)
    ref_length = len(ref_seq)
    feature_start = int(feature.location.start)
    feature_end = int(feature.location.end)
    feature_seq = feature.extract(record_sequence)
    og_feature = deepcopy(feature)
    if 'gene' not in og_feature.qualifiers:
        og_feature.qualifiers['gene'] = [ref_gene_name]
    og_feature_start = int(og_feature.location.start)
    og_feature_end = int(og_feature.location.end)

    def coord_align(ref_seq, feature_seq):
        #Probability that the continuous interval used to find good start/stops
        #occurs by chance should be = 1/(ref_seq*10)
        interval = max(ceil((log(len(ref_seq) * 10))/log(4)), 3)

        #The pairwise aligner returns an iterable list of alignments all with the
        #highest score it could possibly produce. Take the first 5 alignments from this list.
        aligner = Align.PairwiseAligner(scoring="blastn", mode = 'global')
        alignment = aligner.align(ref_seq, feature_seq)
        alignment_list = []
        while len(alignment_list) < 5:
            i = len(alignment_list)
            try:
                target = alignment[i].aligned[0]
                query = alignment[i].aligned[1]
                found_low = (target[0][0] == 0) and (abs(target[0][0] - target[0][1])) >= interval
                found_high = (target[-1][1] == ref_length) and (abs(target[-1][0] - target[-1][1])) >= interval
                alignment_list.append([i, len(target), len(query), found_low, found_high])
                continue
            except IndexError:
                break

        #If possible, filter to find the best fitting alignment rather than the first one available.
        if len(alignment_list) > 1:
            max_found = max(map(lambda x: [x[3],x[4]], alignment_list))
            alignment_list = [i for i in alignment_list if sum([i[3], i[4]]) == sum(max_found)]
            min_breaks = min(map(lambda x: [x[1],x[2]], alignment_list))
            alignment_list = [i for i in alignment_list if sum([i[1], i[2]]) == sum(min_breaks)]
            alignment = alignment[alignment_list[0][0]]
        else:
            alignment = alignment[0]

        target = alignment.aligned[0]
        query = alignment.aligned[1]
        score = alignment.score
        padding = False

        found_low = (target[0][0] == 0) and (abs(target[0][0] - target[0][1])) >= interval
        found_high = (target[-1][1] == ref_length) and (abs(target[-1][0] - target[-1][1])) >= interval

        target_low_seq = ref_seq[target[0][0]:target[0][1]]
        target_high_seq = ref_seq[target[-1][0]:target[-1][1]]
        query_low_seq = feature_seq[query[0][0]:query[0][1]]
        query_high_seq = feature_seq[query[-1][0]:query[-1][1]]

        #make sure there aren't too many mismatches causing falsely assigned found_low/high values
        if (target_low_seq[:3] != query_low_seq[:3]):
            found_low = False
        elif found_low and target[0][1] < (len(ref_seq)/3):
            gaps = ident = mismatch = 0
            for i in range(len(target_low_seq)):
                if target_low_seq[i] == '-' or query_low_seq[i] == '-':
                    gaps += 1
                elif target_low_seq[i] == query_low_seq[i]:
                    ident += 1
                else:
                    mismatch += 1
            if (gaps + mismatch) > (len(target_low_seq)/3):
                found_low = False

        if (target_high_seq[-3:] != query_high_seq[-3:]):
            found_high = False
        elif found_high and abs(target[-1][1] - target[-1][0]) < (len(ref_seq)/3):
            gaps = ident = mismatch = 0
            for i in range(len(target_high_seq)):
                if target_high_seq[i] == '-' or query_high_seq[i] == '-':
                    gaps += 1
                elif target_high_seq[i] == query_high_seq[i]:
                    ident += 1
                else:
                    mismatch += 1
            if (gaps + mismatch) > (len(target_high_seq)/3):
                found_high = False

        if not found_low or not found_high:
            padding = True
        return found_low, found_high, target, query, alignment, padding, score, interval

    def add_padding(feature, target, query, interval):
        #If we're looking to make corrections, add some context to help the aligner
        pad_feature = deepcopy(feature)
        feature_start = feature.location.start
        feature_end = feature.location.end

        #The only time padding matters is when the reference overhangs on the left or right side of
        #the feature. All other scenarios would be where found_low/high = True and extra context is unnecessary.
        left_ref_overhang = int(target[0][0]) > int(query[0][0])
        right_ref_overhang = (int(query[-1][1]) == len(feature)) and (int(target[-1][1]) < len(ref_seq))
        left_pad = right_pad = abs(len(ref_seq) - len(og_feature.extract(record_sequence))) + 2*interval

        if left_ref_overhang:
            left_pad = int(target[0][0] - query[0][0])
            left_pad = left_pad + 2*interval
        if right_ref_overhang:
            right_pad = len(alignment[0]) - int(query[-1][1])
            right_pad = right_pad + 2*interval

        if (fix_start and feature.strand == 1) or (fix_stop and feature.strand == -1):
            if (fix_stop and feature.strand == -1):
                left_pad = right_pad
            padded_feature_start = max(0, (feature_start - left_pad))
            if (not found_low and feature.strand == 1) or (not found_high and feature.strand == -1):
                feature_start = padded_feature_start

        if (fix_stop and feature.strand == 1) or (fix_start and feature.strand == -1):
            if (fix_start and feature.strand == -1):
                right_pad = left_pad
            padded_feature_end = min(len(record_sequence), feature_end + right_pad)
            if (not found_low and feature.strand == -1) or (not found_high and feature.strand == 1):
                feature_end = padded_feature_end

        pad_feature.location = FeatureLocation(
            feature_start,
            feature_end,
            strand=feature.strand
        )
        return pad_feature

    #First alignment
    found_low, found_high, target, query, alignment, padding, first_score, interval = coord_align(ref_seq, feature_seq)
    corrected_feature = deepcopy(feature)
    corrected_feature_start = corrected_feature.location.start
    corrected_feature_end = corrected_feature.location.end
    if padding:
        #Align again after adding padding to the feature sequence if warranted
        pad_feature = add_padding(feature, target, query, interval)
        pad_feature_seq = pad_feature.extract(record_sequence)

        pad_found_low, pad_found_high, pad_target, pad_query, pad_alignment, padding, second_score, second_interval = coord_align(ref_seq, pad_feature_seq)
        if found_low:
            pad_found_low = False
        if found_high:
            pad_found_high = False
    else:
        second_score = -1
        pad_found_high = False
        pad_found_low = False

    for i in range(2):
        if feature.strand == 1:
            good_stop = found_high
            if pad_found_high and fix_stop:
                corrected_feature_end = (pad_feature.location.start + pad_query[-1][1])
                if (second_score > first_score):
                    feature_end = corrected_feature_end
                    good_stop = True
            elif found_high:
                corrected_feature_end = feature_start + query[-1][1]
                if fix_stop:
                    feature_end = corrected_feature_end

            good_start = found_low
            if pad_found_low and fix_start:
                corrected_feature_start = (pad_feature.location.start + (pad_query[0][0]))
                if (second_score > first_score):
                    feature_start = corrected_feature_start
                    good_start = True
            elif found_low:
                corrected_feature_start = feature_start + query[0][0]
                if fix_start:
                    feature_start = corrected_feature_start

        elif feature.strand == -1:
            good_stop = found_high
            if pad_found_high and fix_stop:
                corrected_feature_start = (pad_feature.location.end - pad_query[-1][1])
                if (second_score > first_score):
                    feature_start = corrected_feature_start
                    good_stop = True
            elif found_high:
                corrected_feature_start = feature_end - query[-1][1]
                if fix_stop:
                    feature_start = corrected_feature_start

            good_start = found_low
            if pad_found_low and fix_start:
                corrected_feature_end = (pad_feature.location.end - pad_query[0][0])
                if (second_score > first_score):
                    feature_end = corrected_feature_end
                    good_start = True
            elif found_low:
                corrected_feature_end = feature_end - query[0][0]
                if fix_start:
                    feature_end = corrected_feature_end

        feature.location = FeatureLocation(
            int(feature_start),
            int(feature_end),
            strand=feature.strand
        )
        feature_seq = feature.extract(record_sequence)

        if i == 1:
            continue
        elif any([pad_found_low, pad_found_high]) and any([fix_start, fix_stop]) and (first_score > second_score):
            corrected_feature.location = FeatureLocation(
                int(corrected_feature_start),
                int(corrected_feature_end),
                strand=feature.strand
            )
            corrected_feature_seq = corrected_feature.extract(record_sequence)
            cor_low, cor_high, cor_target, cor_query, cor_alignment, cor_padding, third_score, third_interval = coord_align(ref_seq, corrected_feature_seq)

            if (third_score >= first_score):
                second_score = third_score + 1
                continue
            else:
                break
        else:
            break

    #For cases where we DON'T want to fix coords, good_start/stop should only be true if positionally
    #identical to the start/stop in the reference. This is opposed to the normal definition of
    #good_start/stop as CONTAINING a reference corresponding start/stop site.
    if good_start and not fix_start:
        if ((feature.strand == 1 and feature.location.start != corrected_feature_start) or
            (feature.strand == -1 and feature.location.end != corrected_feature_end)):
            good_start = False
    if good_stop and not fix_stop:
        if ((feature.strand == 1 and feature.location.end != corrected_feature_end) or
            (feature.strand == -1 and feature.location.start != corrected_feature_start)):
            good_stop = False

    if og_feature.location != feature.location:
        feature.qualifiers['translation'] = [
            str(translate(
                feature.extract(record_sequence),
                table=genetic_code,
                to_stop=True,
            ))
        ]
        designator.append_qualifier(
            feature.qualifiers, 'inference',
            "COORDINATES:alignment:Hybran"
        )
        corrected_orf_report.append([og_feature, deepcopy(feature)])
    return good_start, good_stop

def pseudoscan(feature, ref_feature, seq_ident, seq_covg, attempt_rescue=False, blast_hit_dict=None
):
    """
    Determine whether a feature should have the /pseudo qualifier.

    :param feature: SeqFeature object of the one to test for pseudo
    :param ref_feature: SeqFeature object of the reference feature to compare to
    :param seq_ident:
    :param seq_covg:
    :param attempt_rescue: Boolean whether to attempt coordinate correction (feature may still be pseudo after correction)
    :param blast_hit_dict:
        dictionary of blast scores for the reference comparison.
        This will be updated if the coordinates are corrected.
    """
    # Check if the reference annotation was pseudo (some branches of the upcoming conditional use this information)
    # seqret moves the reference pseudo tags to a note field beginning with *pseudo or *pseudogene
    # Save the information and drop the note.
    pseudo_note = False
    if 'note' in feature.qualifiers:
        pseudo_note = [_ for _ in feature.qualifiers['note'] if _.startswith("*pseudo")]
        if pseudo_note:
            feature.qualifiers['note'].remove(pseudo_note[0])

    ref_was_pseudo = pseudo_note or designator.is_pseudo(feature.qualifiers)
    og_broken_stop, og_stop_note = is_broken_stop(feature)
    divisible_by_three = lambda  _: len(_.location) % 3 == 0
    og_feature = deepcopy(feature)
    if ref_was_pseudo:
        good_start, good_stop = coord_check(feature, ref_feature)
        coords_ok = [good_start, good_stop]
        if (not all(coords_ok)) and (divisible_by_three(og_feature) and not divisible_by_three(ref_feature)) and not og_broken_stop:
            feature.qualifiers.pop('pseudo', None)
            feature.qualifiers.pop('pseudogene', None)
            is_pseudo = False
        else:
            if og_stop_note:
                designator.append_qualifier(feature.qualifiers, 'note', og_stop_note)
            is_pseudo = True
            feature.qualifiers['pseudo'] = ['']

    else:
        fix_start = False
        fix_stop = False
        confirmed_feature = False
        og_blast_defined = False
        while not confirmed_feature:
            good_start, good_stop = coord_check(feature, ref_feature, fix_start, fix_stop)
            coords_ok = [good_start, good_stop]
            ref_seq = translate(ref_feature.extract(ref_sequence), table=genetic_code, to_stop=True)
            feature_seq = translate(feature.extract(record_sequence), table=genetic_code, to_stop=True)
            #ref_match with 'thresholds enforced'
            top_hit, low_covg, blast_stats = BLAST.reference_match(
                query=SeqRecord(feature_seq),
                subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                seq_ident=seq_ident,
                seq_covg=seq_covg,
            )
            blast_ok = top_hit and not low_covg
            if og_blast_defined:
                lost_match = (og_blast_ok and not blast_ok)

            broken_stop, stop_note = is_broken_stop(feature)

            if not og_blast_defined:
                if all(coords_ok):
                    confirmed_feature = True
                elif attempt_rescue:
                    fix_start = True
                    # If a gene has lost its stop codon and extended (gene longer than reference),
                    # it will have a "bad_stop" but we will want to keep that bad stop in order to
                    # capture as much information as we can (how far along the gene extended).
                    # However, if it has an early stop (gene is shorter than reference), we want to find out
                    # where it was "supposed" to stop (according to reference) and set the stop position there
                    # to capture as much information as we can.
                    if len(feature.location) < len(ref_feature.location):
                        fix_stop = True
            else:
                # Revert coordinate correction if the blastp hit is lost.
                # (usually due to the reference-corresponding start being affected by an early frameshift)
                if lost_match:
                    feature.location = og_feature.location
                    feature.qualifiers = og_feature.qualifiers
                    broken_stop, stop_note = og_broken_stop, stop_note
                    blast_stats = og_blast_stats
                confirmed_feature = True


            if confirmed_feature:
                if blast_hit_dict:
                    blast_hit_dict.update(blast_stats[ref_feature.qualifiers['gene'][0]])
                if ((all(coords_ok) and divisible_by_three(feature)) or blast_ok) and not broken_stop:
                    is_pseudo = False
                else:
                    is_pseudo = True
                    feature.qualifiers['pseudo']=['']
                    if stop_note:
                        designator.append_qualifier(feature.qualifiers, 'note', stop_note)
                break
            else:
                og_blast_ok, og_blast_stats = (blast_ok, blast_stats)
                og_blast_defined = True

    return is_pseudo

def isolate_valid_ratt_annotations(feature_list, reference_locus_list, seq_ident, seq_covg, ratt_enforce_thresholds,
    nproc=1,
):
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
    non_cds_features = []
    ratt_blast_results = {}
    rejects = []
    valid_features = []

    if ratt_enforce_thresholds:
        ratt_seq_ident = seq_ident
        ratt_seq_covg = seq_covg
    else:
        ratt_seq_ident = ratt_seq_covg = 0

    def refcheck(cds_feature, ratt_seq_ident=ratt_seq_ident, ratt_seq_covg=ratt_seq_covg, record_sequence=record_sequence):
        valid = False
        remark = ''
        blast_stats = {}
        ref_feature = ref_annotation[cds_feature.qualifiers['gene'][0]]
        ref_seq = ref_feature.qualifiers['translation'][0]
        feature_sequence = translate(cds_feature.extract(record_sequence), table=genetic_code, to_stop=True)
        if len(feature_sequence) == 0:
            remark = 'length of AA sequence is 0'
        else:
            top_hit, low_covg, blast_stats = BLAST.reference_match(
                query=SeqRecord(feature_sequence),
                subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                seq_ident=ratt_seq_ident,
                seq_covg=ratt_seq_covg,
            )

            if top_hit:
                valid = True
            else:
                remark = 'No blastp hit to corresponding reference CDS at specified thresholds.'
        return valid, low_covg, blast_stats, remark

    for feature in feature_list:
        if feature.type != 'CDS':
            non_cds_features.append(feature)
            continue

        compound_interval = isinstance(feature.location,Bio.SeqFeature.CompoundLocation)
        # Identify features with 'joins'
        if compound_interval and 'ribosomal_slippage' not in feature.qualifiers:
            #Need to initialize the feature without the compound location attribute.
            #The earliest start and the latest end of the joined feature will be bridged together
            feature_start = feature.location.start
            feature_end = feature.location.end
            feature_strand = feature.strand
            feature.location = (FeatureLocation(feature_start, feature_end, strand=feature_strand))
            #Check if feature has an internal stop codon.
            #
            # If it doesn't, we will assign pseudo and accept it.

            # The gff conversion of a gbk entry with joins is not meaningful,
            # and causes some problems, as the entire sequence gets labeled
            # "biological region" and two basically empty CDS records are created.

            broken_stop, stop_note = is_broken_stop(feature)
            if broken_stop:
                good_start, good_stop = coord_check(
                    feature,
                    ref_annotation[feature.qualifiers['gene'][0]],
                    fix_start=True,
                    fix_stop=True,
                )
                if not good_stop:
                    rejects.append((feature, "RATT-introduced compound interval did not include reference " +
                                    "stop position."))
                    continue
                broken_stop, stop_note = is_broken_stop(feature)
                if broken_stop:
                    feature.qualifiers['pseudo']=['']
                    designator.append_qualifier(feature.qualifiers, 'note', stop_note)
                    valid_features.append(feature)
                    continue
            unbroken_cds.append(feature)
            continue

        feature_is_pseudo = pseudoscan(
            feature,
            ref_annotation[feature.qualifiers['gene'][0]],
            seq_ident,
            seq_covg,
            attempt_rescue=True
        )

        if feature_is_pseudo:
            valid_features.append(feature)
        else:
            unbroken_cds.append(feature)

    logger.debug("Valid CDSs before checking coverage: " + str(len(unbroken_cds) + len(valid_features)))
    logger.debug(f"Checking similarity to reference CDSs using {nproc} process(es)")

    with multiprocessing.Pool(processes=nproc) as pool:
         results = pool.map(
            refcheck,
            unbroken_cds,
        )
    for i in range(len(unbroken_cds)):
        valid, feature_sequence, blast_stats, rejection_note = results[i]
        cds_feature = unbroken_cds[i]
        if valid:
            ratt_blast_results.update(blast_stats)
            cds_feature.qualifiers['translation'] = [str(feature_sequence)]
            valid_features.append(cds_feature)
        else:
            rejects.append((cds_feature, rejection_note))
    logger.debug("Valid CDSs after checking coverage: " + str(len(valid_features)))
    return valid_features, ratt_blast_results, rejects

def find_inframe_overlaps(ratt_features, abinit_features_dictionary):
    """
    Identify the relevant ab initio features.
    Mark for removal those identically annotated by RATT and mark for conflict resolution those with different coordinates, but overlapping in frame.
    :param ratt_features: List of features from RATT (list of SeqFeature objects)
    :param abinit_features_dictionary: sorted dictionary of ab initio features ordered by the genomic position (i.e.,
    feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the corresponding
    SeqFeature object)
    :returns:
        - abinit_features_not_in_ratt (:py:class:`dict`) -
            sorted dictionary of ab initio features that are not annotated by RATT ordered by the genomic position
            i.e. feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the
             corresponding SeqFeature object
        - ratt_overlapping_genes (:py:class:`dict`) -
            dictionary with ab initio feature location triples as keys and a list of RATT feature location triples
            as values. Conflicts here are annotations that overlap in-frame.
        - abinit_rejects (:py:class:`list`) -
            list of tuples of the form (abinit_feature, remark) where abinit_feature is a SeqFeature of a rejected
            ab initio annotation and remark is the rationale as a string.
    """

    abinit_features_not_in_ratt = abinit_features_dictionary.copy()
    ratt_overlapping_genes = collections.defaultdict(list)
    abinit_rejects = []
    for abinit_feature_position in abinit_features_dictionary.keys():
        abinit_feature = abinit_features_dictionary[abinit_feature_position]
        (abinit_start, abinit_end, abinit_strand) = abinit_feature_position
        abinit_duplicate_removed = False
        for ratt_feature in ratt_features:
            if ratt_feature.type != 'CDS':
                continue
            ratt_start = int(ratt_feature.location.start)
            ratt_end = int(ratt_feature.location.end)
            ratt_strand = ratt_feature.location.strand
            if ratt_start > abinit_end and ratt_end > abinit_end:
                break
            elif (overlap_inframe(ratt_feature.location, abinit_feature.location)):
                if len(ratt_feature.location) == len(abinit_feature.location):
                    abinit_rejects.append((abinit_features_not_in_ratt.pop((abinit_start, abinit_end, abinit_strand), None),
                                           'duplicate of ' + ratt_feature.qualifiers['locus_tag'][0]))
                    abinit_duplicate_removed = True
                else:
                    ratt_overlapping_genes[abinit_feature_position].append((ratt_start, ratt_end, ratt_strand))
        if abinit_duplicate_removed:
            ratt_overlapping_genes.pop(abinit_feature_position, None)
    return abinit_features_not_in_ratt, ratt_overlapping_genes, abinit_rejects


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
        mystart = feature.location.start
        myend = feature.location.end
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


def populate_gaps(
        abinit_features,
        intergenic_positions,
        ratt_pre_intergene,
        ratt_post_intergene,
):
    """
    :param abinit_features: dict as from generate_feature_dictionary of ab initio features to try to incorporate
    :param intergenic_positions: (sorted) output of get_interregions
    """

    logger = logging.getLogger('PopulateGaps')
    logger.debug('Merging RATT and ab initio annotations')
    abinit_features = deepcopy(abinit_features)
    abinit_keepers = []
    abinit_conflicts = collections.defaultdict(list)
    prev_unannotated_region_end = {'+':1, '-':1}
    for j in intergenic_positions:
        # Variable definitions
        (ratt_unannotated_region_start, ratt_unannotated_region_end, ratt_strand) = j
        next_abinit_features = deepcopy(abinit_features)
        for feature_position in abinit_features.keys():
            prokka_feature = abinit_features[feature_position]
            (prokka_feature_start, prokka_feature_end, prokka_strand) = feature_position
            ratt_unannotated_region_range = range(ratt_unannotated_region_start,
                                                  ratt_unannotated_region_end + 1)
            if((prokka_strand == -1 and ratt_strand == '-')
               or (prokka_strand == 1 and ratt_strand == '+')):
                # If Prokka feature is location before the end of the previous intergenic region, pop the key from the
                # dictionary and continue loop
                if((prokka_feature_start < prev_unannotated_region_end[ratt_strand]
                    and prokka_feature_end <= prev_unannotated_region_end[ratt_strand])
                   or prokka_feature.type == 'source' # The source feature is accounted for by RATT
                ):
                    next_abinit_features.pop(feature_position, None)
                    continue
                # Else if the prokka feature location is after the end of the intergenic region, break out of
                # the inner loop
                elif(prokka_feature_start > ratt_unannotated_region_end
                     and prokka_feature_end > ratt_unannotated_region_end):
                    break
                # If the ab initio feature is contained in the unannotated range
                elif(prokka_feature_start in ratt_unannotated_region_range
                     and prokka_feature_end in ratt_unannotated_region_range):
                    abinit_keepers.append(prokka_feature)
                # If the Prokka feature overlaps with two RATT features
                elif(prokka_feature_start < ratt_unannotated_region_start
                     and prokka_feature_end > ratt_unannotated_region_end):
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
                    abinit_conflicts[feature_position] += [ratt_overlapping_feature_1_loc, ratt_overlapping_feature_2_loc]
                # If the Prokka feature overlaps with one RATT feature
                else:
                    does_not_overlap = False
                    # the abinit feature might overlap a ratt feature without crossing into the unannotated region.
                    # such cases are more likely to be conflicts, so we should identify them and get them resolved.
                    if (prokka_feature_start < ratt_unannotated_region_start) and \
                            (prokka_feature_end < ratt_unannotated_region_end):
                        ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                       prokka_strand)]
                    elif (prokka_feature_start in ratt_unannotated_region_range) and \
                            prokka_feature_end > ratt_unannotated_region_end:
                        ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                        prokka_strand)]
                    else:
                        does_not_overlap = True
                    if not does_not_overlap:
                        ratt_overlapping_feature_loc = (int(ratt_overlapping_feature.location.start),
                                                        int(ratt_overlapping_feature.location.end),
                                                        int(ratt_overlapping_feature.location.strand))
                        abinit_conflicts[feature_position].append(ratt_overlapping_feature_loc)
        prev_unannotated_region_end[ratt_strand] = ratt_unannotated_region_end
        abinit_features = next_abinit_features

    return abinit_keepers, abinit_conflicts

def check_inclusion_criteria(
        ratt_annotation,
        abinit_annotation,
        reference_gene_locus_dict,
        reference_locus_gene_dict,
        abinit_blast_results,
        ratt_blast_results,
):
    """
    This function compares RATT and Prokka annotations and resolves conflicting annotations.
    Either one feature or both will be accepted.

    :param embl_file:
    :param ratt_annotation:
    :param abinit_annotation:
    :param ratt_gene_location:
    :returns:
        - include_abinit (:py:class:`bool`) - whether the ab initio annotation should be kept
        - include_ratt (:py:class:`bool`) - whether the RATT annotation should be kept
        - remark (:py:class:`str`) - explanation for why the rejected annotation, if any, was not included
    """
    logger = logging.getLogger('CheckInclusionCriteria')
    include_ratt = True
    reject_abinit = False
    remark = ''
    loc_key = (int(abinit_annotation.location.start), int(abinit_annotation.location.end), int(abinit_annotation.location.strand))
    if(abinit_annotation.type != 'CDS'
       or 'gene' not in abinit_annotation.qualifiers.keys()
    ):
        pass
    elif(abinit_annotation.qualifiers['gene'][0] in reference_locus_gene_dict.keys()
         or abinit_annotation.qualifiers['gene'][0] in reference_gene_locus_dict.keys()
    ):
        try:
            locus_tag_list = reference_gene_locus_dict[abinit_annotation.qualifiers['gene'][0]]
        except KeyError:
            locus_tag_list = abinit_annotation.qualifiers['gene']
        locus_tag = ratt_annotation.qualifiers['locus_tag'][0]
        blast_stats = abinit_blast_results[abinit_annotation.qualifiers['locus_tag'][0]]
        same_gene_name = locus_tag in locus_tag_list
        if (not same_gene_name
            and overlap_inframe(abinit_annotation.location, ratt_annotation.location)
            and not designator.is_pseudo(abinit_annotation.qualifiers)
            and not designator.is_pseudo(ratt_annotation.qualifiers)
            ):
            ratt_coord_status = coord_check(ratt_annotation, ref_annotation[ratt_annotation.qualifiers['gene'][0]])
            (ratt_start_ok, ratt_stop_ok) = ratt_coord_status
            ratt_coord_score = sum([int(_) for _ in ratt_coord_status])

            abinit_coord_status = coord_check(abinit_annotation, ref_annotation[abinit_annotation.qualifiers['gene'][0]])
            (abinit_start_ok, abinit_stop_ok) = abinit_coord_status
            abinit_coord_score = sum([int(_) for _ in ratt_coord_status])

            # Both annotations being intact according to their respective reference names
            # suggests that the reference genes are highly similar.
            # RATT's assignment is furthermore based on synteny, so it wins out
            if ((all(ratt_coord_status) and all(abinit_coord_status))
                or ratt_coord_status == abinit_coord_status
                ):
                reject_abinit = True
                remark = f"Equally valid call, but conflicting name with RATT annotation {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)}; RATT favored due to synteny."
            elif (all(ratt_coord_status)
                  or ratt_coord_score > abinit_coord_score
                  or (ratt_stop_ok and not abinit_stop_ok)
                  ):
                reject_abinit = True
                remark = f"RATT annotation {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)} more accurately named and delineated"
            # This is the only other possibility:
            #elif (all(abinit_coord_status)
            #      or abinit_coord_score > ratt_coord_score
            #      or (abinit_stop_ok and not ratt_stop_ok)
            #      ):
            else:
                include_ratt = False
                remark = f"ab initio annotation {extractor.get_ltag(abinit_annotation)}:{extractor.get_gene(abinit_annotation)} more accurately named and delineated"

        elif same_gene_name:
            #Always take the non-pseudo annotation if possible
            if designator.is_pseudo(ratt_annotation.qualifiers) and (not designator.is_pseudo(abinit_annotation.qualifiers)):
                include_abinit = True
                include_ratt = False
                remark = "Non-pseudo ab initio annotation takes precedence."
                return include_abinit, include_ratt, remark
            elif (not designator.is_pseudo(ratt_annotation.qualifiers)) and designator.is_pseudo(abinit_annotation.qualifiers):
                include_abinit = False
                include_ratt = True
                remark = "Non-pseudo ratt annotation takes precedence."
                return include_abinit, include_ratt, remark

            if (locus_tag not in ratt_blast_results.keys()
            ):
                blast_stats = BLAST.reference_match(
                    query=SeqRecord(abinit_annotation.extract(record_sequence)),
                    subject=ref_fna_dict[locus_tag],
                    seq_ident=0,
                    seq_covg=0,
                    blast_type="n"
                )[2]
                blast_stats = blast_stats[locus_tag]
                ratt_blast_results = BLAST.reference_match(
                    query=SeqRecord(ratt_annotation.extract(record_sequence)),
                    subject=ref_fna_dict[locus_tag],
                    seq_ident=0,
                    seq_covg=0,
                    blast_type="n"
                )[2]
            ratt_start = int(ratt_annotation.location.start)
            ratt_stop = int(ratt_annotation.location.end)
            ratt_strand = int(ratt_annotation.location.strand)
            prom_mutation = True
            ratt_prom_seq = upstream_context(ratt_annotation.location, record_sequence)
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
                    prom_mutation = False
                    reject_abinit = True
                    remark = "start position of RATT's " + locus_tag + " corresponds to the reference annotation's"
            if prom_mutation is True:
                ratt_coverage_measure = abs(int(ratt_blast_results[locus_tag]['scov']) -
                                            int(ratt_blast_results[locus_tag]['qcov']))
                prokka_coverage_measure = abs(int(blast_stats['scov']) - int(blast_stats['qcov']))
                if ratt_coverage_measure < prokka_coverage_measure:
                    reject_abinit = True
                    remark = 'RATT annotation for ' + locus_tag + ' has better alignment coverage with the reference'
                elif prokka_coverage_measure < ratt_coverage_measure:
                    logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                    remark = (
                        'Ab initio feature '
                        + abinit_annotation.qualifiers['locus_tag'][0]
                        + ' has better alignment coverage with the reference.'
                    )
                    reject_abinit = False
                    include_ratt = False
                elif ratt_coverage_measure == prokka_coverage_measure:
                    # Checking for identity if coverage is the same
                    if int(ratt_blast_results[locus_tag]['iden']) > int(blast_stats['iden']):
                        reject_abinit = True
                        remark = 'RATT annotation for ' + locus_tag + 'has higher identity with the reference and the same alignment coverage'
                    elif int(blast_stats['iden']) > int(ratt_blast_results[locus_tag]['iden']):
                        logger.debug('Prokka annotation more accurate than RATT for ' + locus_tag)
                        remark = (
                            'Ab initio feature '
                            + abinit_annotation.qualifiers['locus_tag'][0]
                            + ' has higher identity with the reference and the same alignment coverage.'
                        )
                        reject_abinit = False
                        include_ratt = False
                    else:
                        # If RATT and Prokka annotations are same, choose RATT
                        reject_abinit = True
                        remark = 'identical to RATT for ' + locus_tag
                else:
                    reject_abinit = False

    include_abinit = False
    # Check if feature types are the same. If not add feature to EMBL record
    if ratt_annotation.type != abinit_annotation.type:
        include_abinit = True
    # Check if gene names match and if they don't or if gene names are missing, keep both
    elif 'gene' in ratt_annotation.qualifiers.keys() and 'gene' in abinit_annotation.qualifiers.keys():
        if ratt_annotation.qualifiers['gene'] != abinit_annotation.qualifiers['gene']:
            if not reject_abinit:
                include_abinit = True
        # If gene names are the same and the lengths of the genes are comparable between RATT and Prokka annotation
        # (difference in length of less than/equal to 10 bps), the RATT annotation is preferred
        elif ratt_annotation.qualifiers['gene'] == abinit_annotation.qualifiers['gene'] and \
                        abs(len(abinit_annotation.location)-len(ratt_annotation.location)) > 0:
            if not reject_abinit:
                include_abinit = True
        else:
            remark = ("RATT annotation for "
                      + '|'.join([ratt_annotation.qualifiers['locus_tag'][0],ratt_annotation.qualifiers['gene'][0]])
                      + " preferred based on gene length similarity.")
    # If gene tag is missing and the product is not a hypothetical protein, check to see if the products are the
    # same between RATT and Prokka and if they are and if the lengths of the protein coding genes are comparable
    # preferred
    elif ('gene' not in ratt_annotation.qualifiers.keys() or 'gene' not in abinit_annotation.qualifiers.keys()) and \
            ('product' in ratt_annotation.qualifiers.keys() and 'product' in abinit_annotation.qualifiers.keys()):
        if ratt_annotation.qualifiers['product'][0] == 'hypothetical protein' or \
                        abinit_annotation.qualifiers['product'][0] == 'hypothetical protein':
            if not reject_abinit:
                include_abinit = True
        elif ratt_annotation.qualifiers['product'] == abinit_annotation.qualifiers['product'] and \
                        abs(len(abinit_annotation.location)-len(ratt_annotation.location)) > 0:
            if not reject_abinit:
                include_abinit = True
        else:
            remark = ("RATT annotation for product "
                      + '"' + ratt_annotation.qualifiers['product'][0] + '" '
                      + "(" + '|'.join([ratt_annotation.qualifiers['locus_tag'][0],ratt_annotation.qualifiers['gene'][0]]) + ") "
                      + "preferred based on gene length similarity or mismatching product name.")
    else:
        logger.warning('CORNER CASE in check_inclusion_criteria')
        remark = 'corner case: unhandled by inclusion criteria'
    return include_abinit, include_ratt, remark

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


def run(isolate_id, contigs, annotation_fp, ref_proteins_fasta, ref_gbk_fp, reference_genome, script_directory, seq_ident, seq_covg, ratt_enforce_thresholds,
    nproc=1,
):
    """
    Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file> -m
    <output_merged_genes> from the commandline. The log file output stats about the features that are added to the
    RATT annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/'isolate_id'.gbk and
    'isolate_id'/annomerge/'isolate_id'.log

    :param isolate_id: ID of the isolate (Example: H37Rv, 1-0006, etc.). This is the isolate_id that is used for naming
     Genbank files in Prokka
    :param contigs: list of strings for the contig names
    :param annotation_fp: Filepath where RATT, Prokka, reference and prokka no-reference annotations are located.
    Annomerge assumes that RATT annotations are located in <annotation_fp>/ratt, Prokka reference annotations are
    located in <annotation_fp>/prokka and prokka annotations without reference is located in
    <annotation_fp>/prokka-noreference. Additionally annomerge also assumes that withing prokka and prokka-noreference
    directories, the genbank files are located in <isolate_id>.gbk
    :param ref_proteins_fasta: File path for proteome fasta of reference strain
    :param ref_gbk_fp: File path for annotated GenBank file for reference strain
    :param reference_genome: File path for nucleotide fasta of assembled genome
    :param script_dir: Directory where hybran scripts are located
    :param ratt_enforce_thresholds: boolean - whether to enforce seq_ident/seq_covg for RATT-transferred annotations
    :param nproc: int number of processers available for use
    :return: EMBL record (SeqRecord) of annotated isolate
    """

    hybran_tmp_dir = config.hybran_tmp_dir
    global script_dir
    script_dir = script_directory
    global genetic_code
    genetic_code = config.genetic_code
    logger = logging.getLogger('Annomerge')
    logger.debug('Running Annomerge on ' + isolate_id)
    start_time = time.time()

    global reference_gene_locus_dict
    global reference_locus_gene_dict
    reference_gene_list, reference_locus_list, reference_gene_locus_dict, reference_locus_gene_dict, \
        ref_protein_lengths = load_reference_info(ref_proteins_fasta)
    if annotation_fp.endswith('/'):
        file_path = annotation_fp + isolate_id + '/'
    else:
        file_path = annotation_fp + '/' + isolate_id + '/'
    ratt_file_path = file_path + 'ratt'
    ratt_correction_files = []
    ratt_gbk_files = {}
    try:
        for contig in contigs:
            embl_file = f"{isolate_id.replace('|','_')}.{contig.replace('|','_')}.final.embl"
            gbk = converter.convert_embl_to_gbk(ratt_file_path + '/' + embl_file)
            ratt_gbk_files[contig] = gbk
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
    output_merged_genes = os.path.join(isolate_id, 'annomerge', 'merged_genes.gbk')
    output_genbank = os.path.join(isolate_id, 'annomerge', isolate_id + '.gbk')
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    global ratt_rejects  # some RATT annotations are rejected as a side effect of check_inclusion_criteria
    ratt_rejects = []
    ratt_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'ratt_unused.tsv')
    prokka_rejects = []
    prokka_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'prokka_unused.tsv')
    annomerge_records = []
    corrected_abinit_orf_logfile = os.path.join(isolate_id, 'prokka', 'hybran_coord_corrections.tsv')
    corrected_ratt_orf_logfile = os.path.join(isolate_id, 'ratt', 'hybran_coord_corrections.tsv')
    global corrected_orf_report
    corrected_orf_report = []
    ref_contigs = []
    ref_features = []
    # create a dictionary of reference CDS annotations (needed for liftover to ab initio)
    global ref_annotation
    ref_annotation = {}
    # upstream sequence contexts for reference genes. used in multiple places
    global ref_prom_fp_dict
    global ref_fna_dict
    ref_prom_fp_dict = {}
    ref_fna_dict = {}
    for ref_record in SeqIO.parse(ref_gbk_fp, 'genbank'):
        ref_contigs.append(ref_record.seq)
        ref_features.append(ref_record.features)
        prom, fna = get_nuc_seq_for_gene(ref_record.features,ref_record.seq)
        ref_prom_fp_dict.update(prom)
        ref_fna_dict.update(fna)
        for feature in ref_record.features:
            if feature.type != "CDS":
                continue
            # if reference paralogs have been collapsed, the last occurrence in the genome
            # will prevail.
            ref_annotation[feature.qualifiers['gene'][0]] = feature

    output_isolate_recs = []

    for i, contig in enumerate(contigs):
        seqname = '.'.join([isolate_id, contig])
        global ref_sequence
        ref_sequence = ref_contigs[i]
        ratt_contig_record = SeqIO.read(ratt_gbk_files[contig], 'genbank')
        global record_sequence
        record_sequence = ratt_contig_record.seq
        prokka_contig_record = prokka_records[i]

        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        # When prokka assigns the same gene name to multiple orfs, it appends _1, _2, ... to make the names unique.
        # That causes issues for us because we expect all copies of a gene to have the same name.
        for f in prokka_contig_features:
            if 'gene' in f.qualifiers.keys():
                f.qualifiers['gene'][0] = re.sub(r"_\d+$","",f.qualifiers['gene'][0])

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
            # maybe RATT should be adding this inference tag itself
            if 'inference' not in feature.qualifiers:
                feature.qualifiers['inference'] = ["alignment:RATT"]
            else:
                feature.qualifiers['inference'].append("alignment:RATT")
            if feature.type not in [
                    'CDS',
                    'gene',
                    'rRNA', # these aren't reliably transferred by RATT
                    'tRNA',
            ]:
                ratt_contig_non_cds.append(feature)
        logger.debug(f'{seqname}: {len(ratt_contig_non_cds)} non-CDS elements')

        ratt_contig_features, ratt_blast_results, invalid_ratt_features = \
            isolate_valid_ratt_annotations(feature_list=ratt_contig_features,
                                           reference_locus_list=reference_locus_list,
                                           seq_ident=seq_ident,
                                           seq_covg=seq_covg,
                                           ratt_enforce_thresholds=ratt_enforce_thresholds,
                                           nproc=nproc,
            )

        ratt_rejects += invalid_ratt_features
        ratt_contig_features = remove_duplicate_cds(ratt_contig_features)

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
            if corner_cases:
                logger.debug('MERGED GENES: Corner cases')
                for strand in corner_cases_explicit.keys():
                    if len(corner_cases_explicit[strand]) > 0:
                        logger.debug(seqname + ' '.join(corner_cases_explicit[strand]))
        else:
            merged_features = []
        merged_features += identify_conjoined_genes(ratt_contig_features)
        if merged_features and i == 0:
            merged_features_record = prokka_contig_record[:]
            merged_features_record.features = merged_features
            SeqIO.write(merged_features_record, output_merged_genes, 'genbank')
        global ratt_contig_features_dict
        ratt_contig_features_dict = generate_feature_dictionary(ratt_contig_features)
        if len(ratt_contig_features) == 0:
            logger.warning(f"NO RATT ANNOTATION FOR {seqname}")
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
            logger.warning(f"NO AB INITIO ANNOTATION FOR {seqname}")
            prokka_contig_features = ratt_contig_features
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            logger.warning(f'{seqname}: no ab initio annotations to add')
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the
            # features from RATT
            add_prokka_contig_record = prokka_contig_record[:]
            add_prokka_contig_record.features = []
            num_feat = 0
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
            for feat in ratt_contig_record_mod.features:
                if feat.type == 'CDS':
                    added_from_ratt.append(feat.qualifiers['locus_tag'][0])
                else:
                    non_cds_ratt.append(feat)


            abinit_features_dict = generate_feature_dictionary(prokka_contig_features)
            logger.info(f"{seqname}: Checking for in-frame overlaps between RATT and ab initio gene annotations")
            unique_abinit_features_pre_coord_correction, inframe_conflicts_pre_coord_correction, abinit_duplicates = find_inframe_overlaps(
                ratt_contig_features,
                abinit_features_dict,
            )
            prokka_rejects += abinit_duplicates
            logger.debug(f"{seqname}: {len(abinit_duplicates)} ab initio ORFs identical to RATT's")
            logger.debug(f"{seqname}: {len(inframe_conflicts_pre_coord_correction)} ab initio ORFs conflicting in-frame with RATT's")
            logger.debug(f'{seqname}: {len(unique_abinit_features_pre_coord_correction)} total ab initio ORFs remain in consideration')


            logger.debug(f'{seqname}: Checking remaining ab initio CDS annotations for matches to reference using {nproc} process(es)')
            #
            # can contain results for hits to multiple reference genes
            abinit_blast_results_complete = {}
            # only contains results for the accepted reference gene hit
            abinit_blast_results = {}
            refmatch = functools.partial(
                BLAST.reference_match,
                subject=ref_proteins_fasta,
                seq_ident=seq_ident,
                seq_covg=seq_covg,
                identify=lambda _:_.split(':')[1],
            )
            prokka_contig_cdss = [f for f in unique_abinit_features_pre_coord_correction.values() if f.type == 'CDS']
            with multiprocessing.Pool(processes=nproc) as pool:
                blast_package = pool.map(
                    refmatch,
                    [SeqRecord(Seq(f.qualifiers['translation'][0])) for f in prokka_contig_cdss],
                )
            n_coords_corrected = 0
            n_bad_starts = 0
            for j in range(len(prokka_contig_cdss)):
                ref_gene, low_covg, blast_hits = blast_package[j]
                feature = prokka_contig_cdss[j]

                if ref_gene:
                    feature.qualifiers['gene'] = [ref_gene]
                    og_feature_location = deepcopy(feature.location)
                    feature_is_pseudo = pseudoscan(
                        feature,
                        ref_annotation[ref_gene],
                        seq_ident,
                        seq_covg,
                        attempt_rescue=True,
                        blast_hit_dict=blast_hits[ref_gene]
                    )

                    if (og_feature_location != feature.location):
                        n_coords_corrected += 1

                        # for logging purposes
                        if feature_is_pseudo:
                            corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']
                            corrected_orf_report[-1][1].qualifiers['pseudo'] = ['']
                        else:
                            corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']

                    liftover_annotation(
                        feature,
                        ref_annotation[ref_gene],
                        feature_is_pseudo,
                        inference=':'.join([
                            f"similar to AA sequence",
                            os.path.basename(os.path.splitext(ref_gbk_fp)[0]),
                            ref_annotation[ref_gene].qualifiers['locus_tag'][0],
                            ref_gene,
                            "blastp",
                        ])
                    )
                    abinit_blast_results[feature.qualifiers['locus_tag'][0]] = blast_hits[ref_gene]
                # Don't keep gene name assignments from Prokka. They can sometimes be based on
                # poor sequence similarity and partial matches (despite its --coverage option).
                # Keeping them is risky for propagation of the name during clustering.
                elif 'gene' in feature.qualifiers.keys():
                    designator.append_qualifier(
                        feature.qualifiers, 'gene_synonym',
                        feature.qualifiers['gene'][0],
                    )
                    feature.qualifiers.pop('gene', None)
                # We save all the blast hits at this point in case coordinates were corrected and the hits changed
                abinit_blast_results_complete[feature.qualifiers['locus_tag'][0]] = blast_hits

            logger.debug(f'{seqname}: {len(abinit_blast_results.keys())} out of {len(prokka_contig_cdss)} ORFs matched to a reference gene')
            logger.debug(f'{seqname}: Corrected coordinates for {n_coords_corrected} ab initio ORFs')


            # Check for in-frame conflicts/duplicates again since the ab initio gene coordinates changed
            logger.info(f"{seqname}: Checking for fragmented ab initio annotations")
            unique_abinit_features_post_coord_correction_list, dropped_abinit_fragments = process_split_genes(
                unique_abinit_features_pre_coord_correction.values()
            )
            prokka_rejects += dropped_abinit_fragments
            logger.debug(f"{seqname}: {len(dropped_abinit_fragments)} gene fragment pairs merged")
            unique_abinit_features_post_coord_correction = generate_feature_dictionary(
                unique_abinit_features_post_coord_correction_list
            )

            logger.info(f"{seqname}: Checking for new in-frame overlaps with corrected ab initio gene annotations")
            unique_abinit_features, inframe_conflicts, new_abinit_duplicates = find_inframe_overlaps(
                ratt_contig_features,
                unique_abinit_features_post_coord_correction,
            )
            prokka_rejects += new_abinit_duplicates
            logger.debug(f"{seqname}: {len(new_abinit_duplicates)} corrected ab initio ORFs now identical to RATT's")
            logger.debug(f'{seqname}: {len(inframe_conflicts_pre_coord_correction)-len(inframe_conflicts)} corrected ab initio ORFs now conflicting in-frame with RATT')
            logger.debug(f'{seqname}: {len(unique_abinit_features)} total ab initio ORFs remain in consideration')


            intergenic_ratt, intergenic_positions, ratt_pre_intergene, ratt_post_intergene = \
                get_interregions(ratt_contig_record_mod, intergene_length=1)
            sorted_intergenic_positions = sorted(intergenic_positions)
            add_features_from_prokka, overlap_conflicts = populate_gaps(
                abinit_features=unique_abinit_features,
                intergenic_positions=sorted_intergenic_positions,
                ratt_pre_intergene=ratt_pre_intergene,
                ratt_post_intergene=ratt_post_intergene,
            )
            add_prokka_contig_record.features += add_features_from_prokka
            logger.debug(f"{seqname}: {len(add_features_from_prokka)} ab initio ORFs fall squarely into RATT's CDS-free regions.")

            # merge the two dicts of lists. hat tip to https://stackoverflow.com/a/5946322
            abinit_conflicts = collections.defaultdict(list)
            for abinit_conflict_group in (inframe_conflicts, overlap_conflicts):
                for key, value in abinit_conflict_group.items():
                    abinit_conflicts[key] += value

            logger.debug(f"{seqname}: {len(abinit_conflicts.keys())} ab initio CDSs in total overlap RATT CDSs. Resolving...")

            for feature_position in abinit_conflicts.keys():
                abinit_feature = unique_abinit_features[feature_position]
                # Conflict Resolution
                for ratt_conflict_loc in abinit_conflicts[feature_position]:
                    # if the RATT annotation got rejected at some point, its remaining conflicts are moot
                    if ratt_conflict_loc not in ratt_contig_features_dict.keys():
                        include_abinit = True
                        continue
                    ratt_feature = ratt_contig_features_dict[ratt_conflict_loc]
                    # When the ab initio feature is in-frame with a RATT transferred gene
                    # but didn't itself get a reference gene name assigned via direct comparison,
                    # we check whether it has a match to the reference-transferred gene itself
                    # and assign the name if it does. If the name is established, the conflict resolution
                    # is better informed.
                    if(abinit_feature.type == 'CDS'
                       and 'gene' not in abinit_feature.qualifiers.keys()
                       and overlap_inframe(ratt_feature.location, abinit_feature.location)
                       ):
                        if 'translation' not in ratt_contig_features_dict[ratt_conflict_loc].qualifiers.keys():
                            blast_type = "n"
                            query = SeqRecord(abinit_feature.extract(record_sequence))
                            subject = SeqRecord(
                                ratt_contig_features_dict[ratt_conflict_loc].extract(record_sequence),
                                id=ratt_contig_features_dict[ratt_conflict_loc].qualifiers['gene'][0],
                            )
                        else:
                            blast_type = "p"
                            query=SeqRecord(Seq(abinit_feature.qualifiers['translation'][0]))
                            subject=SeqRecord(
                                Seq(ratt_contig_features_dict[ratt_conflict_loc].qualifiers['translation'][0]),
                                id=ratt_contig_features_dict[ratt_conflict_loc].qualifiers['gene'][0],
                            )
                        ref_gene = BLAST.reference_match(
                            query=query,
                            subject=subject,
                            blast_type=blast_type,
                            seq_ident=seq_ident,
                            seq_covg=seq_covg,
                        )[0]
                        if ref_gene:
                            # reference_match's pseudo determination may not be accurate since we are here
                            # comparing annotations on the same genome, so
                            # a full-length match to an already truncated gene should not be considered
                            # a complete reference match
                            abinit_blast_results[abinit_feature.qualifiers['locus_tag'][0]] = \
                                abinit_blast_results_complete[abinit_feature.qualifiers['locus_tag'][0]][ref_gene]
                            pseudo = pseudoscan(
                                feature,
                                ref_feature=ref_annotation[ref_gene],
                                seq_ident=seq_ident,
                                seq_covg=seq_covg,
                                attempt_rescue=True,
                                blast_hit_dict=abinit_blast_results[abinit_feature.qualifiers['locus_tag'][0]]
                            )

                            liftover_annotation(
                                abinit_feature,
                                ref_annotation[ref_gene],
                                pseudo,
                                inference=':'.join([
                                    "similar to AA sequence",
                                    os.path.basename(os.path.splitext(ref_gbk_fp)[0]),
                                    ref_annotation[ref_gene].qualifiers['locus_tag'][0],
                                    ref_gene,
                                    "RATT+blastp",
                                ])
                            )
                    include_abinit, include_ratt, remark = check_inclusion_criteria(
                        ratt_annotation=ratt_feature,
                        abinit_annotation=abinit_feature,
                        reference_gene_locus_dict=reference_gene_locus_dict,
                        reference_locus_gene_dict=reference_locus_gene_dict,
                        abinit_blast_results=abinit_blast_results,
                        ratt_blast_results=ratt_blast_results,
                    )
                    if not include_abinit:
                        prokka_rejects.append((abinit_feature,remark))
                        break
                    elif not include_ratt:
                        ratt_rejects.append((ratt_contig_features_dict.pop(ratt_conflict_loc), remark))
                # Add the abinit feature if it survived all the conflicts
                if include_abinit:
                    add_prokka_contig_record.features.append(abinit_feature)

            if len(merged_features) > 0:
                for feature_1 in merged_features:
                    add_prokka_contig_record.features.append(feature_1)
            for ratt_feature_append in ratt_contig_features_dict.values():
                add_prokka_contig_record.features.append(ratt_feature_append)
            for non_cds in ratt_contig_non_cds:
                add_prokka_contig_record.features.append(non_cds)

            annomerge_records.append(add_prokka_contig_record)

        # Post-processing of genbank file to remove duplicates and rename locus_tag for
        # Prokka annotations
        seqname = '.'.join([isolate_id, contig])
        annomerge_records[i].name = seqname
        # TODO - replace version variable with importlib.version call (and probably url too) in python 3.8+
        annomerge_records[i].annotations['comment'] = "Annotated using hybran " + __version__ + " from https://lpcdrp.gitlab.io/hybran."
        prokka_rec = annomerge_records[i]
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
                if 'Fusion' in feature.qualifiers['note'][0]:
                    if 'gene' in feature.qualifiers.keys():
                        feature.qualifiers.pop('gene')
                    loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                    added_cds[loc_key] = feature
                    prokka_rec.features.append(feature)
            if feature.type == 'CDS':
                loc_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                if loc_key not in added_cds.keys():
                    added_cds[loc_key] = feature
                    prokka_rec.features.append(feature)
                elif added_cds[loc_key].qualifiers['locus_tag'] == feature.qualifiers['locus_tag']:
                    continue
                else:
                    prokka_rec.features.append(feature)
        # Adding translated sequence to RATT annotations
        for feature in prokka_rec.features:
            if feature.type == 'CDS' and 'translation' not in feature.qualifiers.keys():
                feature_sequence = translate(feature.extract(record_sequence), table=genetic_code, to_stop=True)
                feature.qualifiers['translation'] = [feature_sequence]

        logger.debug(f'{seqname}: final feature annotation verification')
        added_ltags = []
        for feature_final in prokka_rec.features:
            if feature_final.type != 'CDS':
                continue
            if 'gene' not in feature_final.qualifiers.keys():
                feature_final.qualifiers['gene'] = feature_final.qualifiers['locus_tag']
            added_ltags.append(feature_final.qualifiers['locus_tag'][0])
        sorted_final = get_ordered_features(prokka_rec.features)
        prokka_rec.features = sorted_final
        output_isolate_recs.append(prokka_rec)
        isolate_features = prokka_rec.features
        output_isolate_recs[i].features = []
        prev_feature_list = []
        num_overlaps = 0
        positions_to_be_resolved = []
        resolve_pairs = []
        for feature in isolate_features:
            if feature.type != 'CDS':
                continue
            if designator.is_raw_ltag(feature.qualifiers['locus_tag'][0]):
                feature.hybran_supplier = 'abinit'
            else:
                feature.hybran_supplier = 'ratt'
            if len(prev_feature_list) == 0:
                prev_feature_list = [int(feature.location.start), int(feature.location.end),
                                     int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                prev_feature = feature
                continue
            if (feature.hybran_supplier != prev_feature.hybran_supplier
                and ((int(feature.location.start) <= prev_feature_list[1]
                      and int(feature.location.start) >= prev_feature_list[0]
                      and (int(feature.location.start) == prev_feature_list[0] or
                           int(feature.location.end) == prev_feature_list[1])
                      )
                     or
                    (int(feature.location.end) <= prev_feature_list[1]
                     and int(feature.location.end) >= prev_feature_list[0]
                     and (int(feature.location.start) == prev_feature_list[0] or
                          int(feature.location.end) == prev_feature_list[1])
                    ))
                ):
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
                output_isolate_recs[i].features.append(feature)
            else:
                position = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
                if position not in positions_to_be_resolved:
                    output_isolate_recs[i].features.append(feature)
        for feat_pair in resolve_pairs:
            if designator.is_raw_ltag(feat_pair[0].qualifiers['locus_tag'][0]):
                prokka_annotation = feat_pair[0]
                ratt_annotation = feat_pair[1]
            else:
                prokka_annotation = feat_pair[1]
                ratt_annotation = feat_pair[0]
            take_abinit, take_ratt, remark = check_inclusion_criteria(
                ratt_annotation=ratt_annotation,
                abinit_annotation=prokka_annotation,
                reference_gene_locus_dict=reference_gene_locus_dict,
                reference_locus_gene_dict=reference_locus_gene_dict,
                abinit_blast_results=abinit_blast_results,
                ratt_blast_results=ratt_blast_results,
            )
            if take_ratt:
                output_isolate_recs[i].features.append(ratt_annotation)
            else:
                ratt_rejects.append((ratt_annotation, remark))
            if take_abinit:
                output_isolate_recs[i].features.append(prokka_annotation)
            else:
                prokka_rejects.append((prokka_annotation, remark))
        ordered_feats = get_ordered_features(output_isolate_recs[i].features)
        # Remove AA translation from pseudos
        for feature in ordered_feats:
            if designator.is_pseudo(feature.qualifiers):
                feature.qualifiers.pop('translation', None)

        output_isolate_recs[i].features = ordered_feats[:]
        final_cdss = [f for f in ordered_feats if f.type == 'CDS']
        logger.debug(f'{seqname}: {len(final_cdss)} CDSs annomerge')

    with open(ratt_rejects_logfile, 'w') as ratt_rejects_log:
        [log_feature_fate(_[0], ratt_rejects_log, _[1]) for _ in ratt_rejects]
    with open(prokka_rejects_logfile, 'w') as prokka_rejects_log:
        # we sometimes get the same gene appearing here multiple times, so collapse
        # the list and merge the remarks.
        prokka_rejects.sort(key=lambda _:_[0].qualifiers['locus_tag'][0])
        last_gene = prokka_rejects[0][0]
        last_remark = ''
        for (gene, remark) in prokka_rejects:
            if gene == last_gene:
                if last_remark:
                    remark = ';'.join([last_remark,remark])
            else:
                log_feature_fate(last_gene, prokka_rejects_log, last_remark)
            last_gene = gene
            last_remark = remark
        log_feature_fate(last_gene, prokka_rejects_log, last_remark)

    with open(corrected_abinit_orf_logfile, 'w') as abinit_corlog, \
         open(corrected_ratt_orf_logfile, 'w') as ratt_corlog:
        header = ['locus_tag', 'gene_name', 'strand', 'og_start', 'og_end', 'new_start', 'new_end', 'fixed_start_codon', 'fixed_stop_codon', 'status']
        print('\t'.join(header), file=abinit_corlog)
        print('\t'.join(header), file=ratt_corlog)
        for (orig_feature, corr_feature) in corrected_orf_report:
            if designator.is_raw_ltag(orig_feature.qualifiers['locus_tag'][0]):
                logfile = abinit_corlog
            else:
                logfile = ratt_corlog
            log_coord_correction(orig_feature,
                                 corr_feature,
                                 logfile,
            )

    SeqIO.write(output_isolate_recs, output_genbank, 'genbank')

    logger.debug('annomerge run time: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes')
