#!/usr/bin/env python2.7

import os
import pickle
import re
import sys
import tempfile
import argparse
import glob
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastpCommandline


GROUPHOME = os.environ['GROUPHOME']


def create_isolate_seq_dict(directory):
    isolate_seq_dict = {}
    for genbank_file in os.listdir(directory):
        if genbank_file.endswith('.gbk'):
            for record in SeqIO.parse(genbank_file, 'genbank'):
                features = record.features
                break
            isolate = genbank_file.split('.')[0]
            isolate_seq_dict[isolate] = {}
            for feature in features:
                if feature.type == 'CDS':
                    seq = feature.qualifiers['translation'][0]
                    isolate_seq_dict[isolate][feature.qualifiers['locus_tag'][0]] = seq
    return isolate_seq_dict


def parse_clustered_proteins(clustered_proteins, annotations):
    underscore_re = re.compile('_[0-9]$')

    # Private Function #
    ##################################################################################################################
    def gff_dict(annotation_dir):
        gff_dictionary = {}
        for gff in os.listdir(annotation_dir):
            isolate_id_ltag = {}
            try:
                with open(annotation_dir + gff, 'r') as gff_file:
                    for line in gff_file:
                        if line.startswith('#'):
                            continue
                        gene = ''
                        column = line.rstrip('\n').split('\t')
                        if len(column) >= 8:
                            info = column[8].split(';')
                            gene_id = ''.join([i.split('=')[1] for i in info if i.startswith('ID=')])
                            locus_tag = ''.join([i.split('=')[1] for i in info if i.startswith('locus_tag=')])
                            gene = ','.join([i.split('=')[1] for i in info if i.startswith('gene=')])
                            if not locus_tag.startswith('Rv') and not locus_tag.startswith('L'):
                                gene = locus_tag
                                locus_tag = ''
                            isolate_id_ltag[gene_id] = (locus_tag, gene)
                    gff_dictionary[isolate] = isolate_id_ltag
            except IOError:
                continue
        return gff_dictionary

    # Different gene names/Rvs in cluster
    def different_genes(gene_cluster_list):
        genes = list(set(sorted([gene[1] for gene in gene_cluster_list])))
        locus_tags = list(set(sorted([locus[0] for locus in gene_cluster_list if locus[0]])))
        return genes, locus_tags

    def find_underscores(gene_cluster_list):
        for tup in gene_cluster_list:
            if underscore_re.search(tup[1]):
                return True

    ##################################################################################################################

    gffs = gff_dict(annotations)
    representative_fasta_list = []
    gene_cluster = OrderedDict()
    different_genes_cluster_w_ltags = {}
    same_genes_cluster_w_ltags = {}
    underscores = {}
    l_tag_only_clusters = {}
    unique_genes_list = []
    with open(clustered_proteins, 'r') as clustered_proteins_file:
        for line in clustered_proteins_file:
            isolates_ids = line.rstrip('\n').split('\t')

            representative_seq_id = isolates_ids[0].split(': ')[1]
            rep_isolate = representative_seq_id.split('-L')[0]
            rep_seq_id = 'L' + representative_seq_id.split('-L')[1]
            representative_ltag_gene_tup = gffs[rep_isolate][rep_seq_id]
            representative = ','.join([rep_isolate, ','.join(representative_ltag_gene_tup)])
            # Getting clusters with only 1 gene that is an L tag or an underscore
            if len(isolates_ids) == 1 and \
                    (representative_ltag_gene_tup[0].startswith('L') or
                     underscore_re.search(representative_ltag_gene_tup[0])):
                unique_genes_list.append([rep_isolate] + list(representative_ltag_gene_tup))

            representative_fasta_list.append([rep_isolate] + list(representative_ltag_gene_tup))
            gene_cluster[representative] = []

            cluster_list = []
            cluster_list_w_isolate = []
            cluster_list.append(representative_ltag_gene_tup)
            cluster_list_w_isolate.append([rep_isolate] + list(representative_ltag_gene_tup))
            for isolate_gene_id in isolates_ids[1:]:
                isolate = isolate_gene_id.split('-L')[0]
                gene_id = 'L' + isolate_gene_id.split('-L')[1]
                cluster_list.append(gffs[isolate][gene_id])
                cluster_list_w_isolate.append([isolate] + list(gffs[isolate][gene_id]))
                gene_cluster[representative].append(','.join([isolate, ','.join(gffs[isolate][gene_id])]))
            number_genes, number_locus_tags = different_genes(cluster_list)
            unique_genes = []
            unique_locus = []
            for data in number_genes:
                if not data.startswith('L'):
                    unique_genes.append(data)
            for data in number_locus_tags:
                if not data.startswith('L'):
                    unique_locus.append(data)
            # All different with L tags
            if (len(unique_genes) > 1 or len(unique_locus) > 1) and \
                    any([gene[0].startswith('L') for gene in cluster_list]):
                different_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            #
            elif (len(unique_genes) == 1 or len(unique_locus) == 1) and \
                    any([gene[0].startswith('L') for gene in cluster_list]):
                same_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            # Getting genes with underscores (partial hits)
            if find_underscores(cluster_list):
                underscores[representative] = cluster_list_w_isolate
            # L tag only clusters
            if all(locus[0].startswith('L') for locus in cluster_list):
                l_tag_only_clusters[representative] = cluster_list_w_isolate
    return different_genes_cluster_w_ltags, same_genes_cluster_w_ltags, underscores, l_tag_only_clusters, \
        unique_genes_list


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


def identify_top_hits(blast_output_file, identity=95, coverage=95):
    """
    :param blast_output_file: Output file from Blastp
    :return: Dictionary of top hits for each unannotated locus that passes the identity and coverage threshold
    """
    all_hits_dict = {}
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
        else:
            continue
    if rv_hit:
        top_hit = get_top_hit(all_hits_dict)
    else:
        top_hit = None
    return top_hit, all_hits_dict


def update_dictionary_ltag_assignments(isolate_id, isolate_ltag, new_gene_name):
    logger = logging.getLogger('GeneNameAssignment')
    if isolate_id not in isolate_update_dictionary.keys():
        isolate_update_dictionary[isolate_id] = {}
        isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
    else:
        isolate_dict_added = isolate_update_dictionary[isolate_id]
        if isolate_ltag in isolate_dict_added.keys():
            logger.warn('ERROR: This locus tag already has an Rv' + '\n')
            logger.warn('Locus Tag:' + '\n')
            logger.warn(isolate_id + '\n')
            logger.warn(isolate_ltag + '\n')
            logger.warn(isolate_update_dictionary[isolate_id][isolate_ltag] + '\n')
            logger.warn('Above locus tag is also annotated as: ' + new_gene_name + '\n')
        else:
            isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
    return


def get_cluster_fasta(rep, cluster_list, isolates_dir):
    rep_isolate_id = rep.split(',')[0]
    rep_locus = rep.split(',')[1]
    rep_gene_name = rep.split(',')[2]
    cluster_fasta_fp = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
    cluster_temp_fasta = cluster_fasta_fp.name
    added_seq = []
    unannotated_genes_list = []
    fasta_tmp = open(cluster_temp_fasta, 'w')
    if not (rep_locus.startswith('L') and rep_gene_name.startswith('L')):
        rep_isolate_fp = isolates_dir + '/' + rep_isolate_id + '.gbk'
        rep_sequence = isolate_sequences[rep_isolate_fp][rep_locus]
        rep_sequence_object = Seq(rep_sequence)
        seq_id = '|'.join([rep_isolate_id, rep_locus, rep_gene_name])
        added_seq.append(seq_id)
        rep_record = SeqRecord(rep_sequence_object, id=seq_id)
        SeqIO.write(rep_record, fasta_tmp, 'fasta')
    else:
        rep_list = [rep_isolate_id, rep_locus, rep_gene_name]
        unannotated_genes_list.append(rep_list)
    for gene in cluster_list:
        gene_isolate_id = gene[0]
        gene_locus = gene[1]
        gene_name = gene[2]
        if gene_locus.startswith('L') and gene_name.startswith('L'):
            unannotated_genes_list.append(gene)
        else:
            gene_id = '|'.join([gene_isolate_id, gene_locus, gene_name])
            if gene_id in added_seq:
                continue
            else:
                gene_isolate_fp = isolates_dir + '/' + gene_isolate_id + '.gbk'
                gene_sequence = isolate_sequences[gene_isolate_fp][gene_locus]
                gene_sequence_object = Seq(gene_sequence)
                added_seq.append(gene_id)
                gene_record = SeqRecord(gene_sequence_object, id=gene_id)
                SeqIO.write(gene_record, fasta_tmp, 'fasta')
    print(len(added_seq))
    return cluster_temp_fasta, unannotated_genes_list


def cluster_annotation_presence(cluster_list):
    # Assumes representative sequence is also present in cluster_list
    annotated = 0
    for gene in cluster_list:
        gene_locus = gene[1]
        gene_name = gene[2]
        if not (gene_locus.startswith('L') and gene_name.startswith('L')):
            annotated += 1
    if annotated == len(cluster_list):
        return True
    else:
        return False


def find_larges_mtb_increment(annotation_directory):
    args = ['grep', 'gene=MTB']
    args.extend(glob.glob(annotation_directory + '/*.gff'))
    stdout = subprocess.Popen(args,
                              stdout=subprocess.PIPE)
    largest_mtb = sorted([i.split('=')[1] for line in stdout.stdout
                          for i in line.rstrip('\n').split(';') if i.startswith('gene=MTB')])[-1]
    if largest_mtb:
        return largest_mtb
    return 0


def ref_seqs(gbk_dir):
    nucleotide_cds = []
    protein_cds = []
    for gbk in os.listdir(gbk_dir):
        if gbk.endswith('.gbk'):
            gbk_filename = gbk.split('.')[0]
            for record in SeqIO.parse(gbk_dir + gbk, 'genbank'):
                if record.features:
                    for feature in record.features:
                        if feature.type == 'CDS':
                            seq = SeqRecord(feature.location.extract(record).seq,
                                            id=feature.qualifiers['gene'][0],
                                            description=gbk_filename)
                            nucleotide_cds.append(seq)
                            seq = SeqRecord(feature.qualifiers['translation'],
                                            id=feature.qualifiers['gene'][0],
                                            description=gbk_filename)
                            protein_cds.append(seq)
    return nucleotide_cds, protein_cds


def blast(subject, stdin_seq):
    blast_command = NcbiblastpCommandline(subject=subject,
                                          outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
    stdout, stderr = blast_command(stdin=stdin_seq)
    return stdout


def find_unannotated_genes(reference_protein_fasta):
    unannotated_seqs = []
    for record in SeqIO.parse(reference_protein_fasta, 'fasta'):
        if record.id.startswith('MTB'):
            unannotated_seqs.append(record)
    with open('unannotated_seqs.fasta', 'w') as out:
        for s in unannotated_seqs:
            SeqIO.write(s, out, 'fasta')
    return 'unannotated_seqs.fasta'


def singleton_clusters(singleton_dict, annotation_dir, reference_fasta, unannotated_fasta, mtb_increment):
    logger = logging.getLogger('SingletonClusters')
    # 3. If a cluster has only one gene (with no gene names), then BLAST representative sequence to H37Rv and if there
    # is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster is annotated
    #  with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes in the
    # cluster.
    logger.info('Number of singleton clusters with single genes: ' + str(len(singleton_dict)) + '\n')
    for single_gene in singleton_dict:
        isolate_id = single_gene[0]
        locus_tag = single_gene[1]
        gene_name = single_gene[2]
        isolate_fp = annotation_dir + '/' + isolate_id + '.gbk'
        if gene_name.startswith('L') and locus_tag.startswith('L'):
            gene_sequence = isolate_sequences[isolate_fp][locus_tag]
            stdout = blast(reference_fasta, gene_sequence)
            top_hit, all_hits = identify_top_hits(stdout)
            assign_mtb = False
            name_to_assign = ''
            if top_hit is None:
                assign_mtb = True
                if os.path.exists(unannotated_fasta):
                    fasta_records = SeqIO.parse(unannotated_fasta, 'fasta')
                    for mtb_record in fasta_records:
                        mtb_id = mtb_record.description.split(' ')[0]
                        mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))
                    stdout_2 = blast(unannotated_fasta, gene_sequence)
                    top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
                    if top_hit_mtb is None:
                        assign_mtb = True
                    else:
                        assign_mtb = False
                        name_to_assign = top_hit_mtb
            else:
                name_to_assign = top_hit
            if assign_mtb:
                mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                mtb_increment = mtb_increment + 1
                sequence_object = Seq(gene_sequence)
                seq_record = SeqRecord(sequence_object, id=mtb_id)
                with open(unannotated_fasta, 'a') as mtb_fasta:
                    SeqIO.write(seq_record, mtb_fasta, 'fasta')
                name_to_assign = mtb_id
                logger.info('Assigned new mtb id' + '\n')
                logger.info(name_to_assign + '\n')
            else:
                logger.info('Assigned existing Rv or MTB' + '\n')
                logger.info(name_to_assign + '\n')
            update_dictionary_ltag_assignments(isolate_id, locus_tag, name_to_assign)
        else:
            continue
    return mtb_increment


def only_ltag_clusters(in_dict, annotation_dir, reference_fasta, unannotated_fasta, mtb_increment):
    logger = logging.getLogger('NewGeneClusters')
    # 2. If a cluster has only candidate novel genes (with no gene names), then BLAST representative sequence to H37Rv
    # and if there is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster
    # is annotated with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes
    # in the cluster.
    logger.info('Number of clusters with only L-tags genes: ' + str(len(in_dict.keys())) +
                '\n')
    for rep_gene in in_dict.keys():
        rep_isolate_id = rep_gene.split(',')[0]
        rep_locus = rep_gene.split(',')[1]
        rep_gene_name = rep_gene.split(',')[2]
        rep_isolate_fp = annotation_dir + '/' + rep_isolate_id + '.gbk'
        rep_sequence = isolate_sequences[rep_isolate_fp][rep_locus]
        # Writing representative amino acid sequence to a temp file to check if all genes in cluster share identity with
        #  this sequence
        rep_fp = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
        rep_temp_fasta = rep_fp.name
        rep_sequence_object = Seq(rep_sequence)
        rep_record = SeqRecord(rep_sequence_object, id='|'.join([rep_isolate_id, rep_locus, rep_gene_name]))
        with open(rep_temp_fasta, 'w') as fasta_tmp:
            SeqIO.write(rep_record, fasta_tmp, 'fasta')
        # Blasting representative sequence with H37Rv
        stdout = blast(reference_fasta, rep_sequence)
        top_hit, all_hits = identify_top_hits(stdout)
        assign_mtb = False
        name_to_assign = ''
        if top_hit is None:
            assign_mtb = True
            if os.path.exists(unannotated_fasta):
                fasta_records = SeqIO.parse(unannotated_fasta, 'fasta')
                for mtb_record in fasta_records:
                    mtb_id = mtb_record.description.split(' ')[0]
                    mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))
                stdout_2 = blast(unannotated_fasta, rep_sequence)
                top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
                if top_hit_mtb is None:
                    assign_mtb = True
                else:
                    assign_mtb = False
                    name_to_assign = top_hit_mtb
        else:
            name_to_assign = top_hit
        if assign_mtb:
            mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
            mtb_increment = mtb_increment + 1
            sequence_object = Seq(rep_sequence)
            seq_record = SeqRecord(sequence_object, id=mtb_id)
            with open(unannotated_fasta, 'a') as mtb_fasta:
                SeqIO.write(seq_record, mtb_fasta, 'fasta')
            name_to_assign = mtb_id
            logger.info('Assigned new mtb id' + '\n')
            logger.info(name_to_assign + '\n')
        else:
            logger.info('Assigned existing Rv or MTB' + '\n')
            logger.info(name_to_assign + '\n')
        for gene_in_cluster in in_dict[rep_gene]:
            update_dictionary_ltag_assignments(gene_in_cluster[0], gene_in_cluster[1], name_to_assign)
        os.unlink(rep_temp_fasta)
    return mtb_increment


def single_gene_clusters(single_gene_dict):
    logger = logging.getLogger('SingleGeneClusters')
    # 1. Handling clusters with only one Rv/gene name in cluster
    # 1. If cluster has candidate novel genes (L_***** or L2_*****) clustered together with a H37Rv annotation, update
    #  all candidate novel genes in this cluster with the H37Rv gene name.
    rep_ltag_keys = []
    logger.info('Number of clusters with single genes: ' + str(len(single_gene_dict.keys())) + '\n')
    for key_sgc in single_gene_dict:
        rep_ltag = False
        key_elements_sgc = key_sgc.split(',')
        if key_elements_sgc[1].startswith('Rv'):
            gene_to_add = key_elements_sgc[1]
        elif key_elements_sgc[2].startswith('Rv'):
            gene_to_add = key_elements_sgc[2]
        elif key_elements_sgc[1].startswith('L') and key_elements_sgc[2].startswith('L'):
            rep_ltag = True
            rep_ltag_keys.append(key_sgc)
        else:
            gene_to_add = key_elements_sgc[2]
        # If representative sequence is an Rv and not an L-tag
        if not rep_ltag:
            genes_in_cluster = single_gene_dict[key_sgc]
            for gene in genes_in_cluster:
                if gene[1].startswith('L') and gene[2].startswith('L'):
                    update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                else:
                    continue
        # If representative sequence is a L-tag
        else:
            genes_in_cluster = single_gene_dict[key_sgc]
            gene_to_add = ''
            for gene in genes_in_cluster:
                if gene[1].startswith('L') and gene[2].startswith('L'):
                    continue
                else:
                    if not gene[1].startswith('L'):
                        gene_to_add = gene[1]
                    else:
                        gene_to_add = gene[2]
                    break
            if len(gene_to_add) == 0:
                logger.info('No Rv found in this cluster' + '\n')
                logger.info(key_sgc + '\n')
                logger.info(genes_in_cluster + '\n')
            else:
                for gene in genes_in_cluster:
                    if gene[1].startswith('L') and gene[2].startswith('L'):
                        update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                    else:
                        continue


def multigene_clusters(in_dict, single_gene_cluster_complete, annotation_dir, unannotated_fasta, mtb_increment):
    logger = logging.getLogger('MultigeneClusters')
    # 4. If a cluster has multiple H37Rv genes and L_tags, BLAST the L_tags to all genes in the cluster that are H37Rv
    # and annotate L_tag with the top hit.
    mgc_output = []

    num_multi = 0
    for gene in in_dict.keys():
        num_multi += 1
        output_line = ''
        unassigned_l_tags = []
        true_multi_cluster = False
        gene_elements = gene.split(',')
        for gene_in_cluster in in_dict[gene]:
            if gene_in_cluster[1] in gene_elements or gene_in_cluster[2] in gene_elements:
                continue
            elif gene_in_cluster[1].startswith('L') and gene_in_cluster[2].startswith('L'):
                unassigned_l_tags.append(gene_in_cluster)
            else:
                true_multi_cluster = True
        if true_multi_cluster:
            output_line = output_line + gene
            for gene_in_cluster in in_dict[gene]:
                gene_str = ','.join(gene_in_cluster)
                output_line = output_line + '\t' + gene_str
            output_line = output_line + '\n'
            mgc_output.append(output_line)
            all_genes_annotated = cluster_annotation_presence(in_dict[gene])
            # If all genes in the cluster are annotated, do nothing. If there is an unannotated gene in the cluster,
            # blast it to all annotated genes in cluster and annotated with top hit
            if not all_genes_annotated:
                fasta_fp_to_blast, genes_to_annotate = get_cluster_fasta(gene, in_dict[gene], annotation_dir)
                for unannotated_gene in genes_to_annotate:
                    unannotated_gene_isolate = unannotated_gene[0]
                    unannotated_gene_locus = unannotated_gene[1]
                    unannotated_isolate_fp = annotation_dir + '/' + unannotated_gene_isolate + '.gbk'
                    unannotated_gene_seq = isolate_sequences[unannotated_isolate_fp][unannotated_gene_locus]
                    stdout = blast(fasta_fp_to_blast, unannotated_gene_seq)
                    top_hit, all_hits = identify_top_hits(stdout)
                    assign_mtb = False
                    name_to_assign = ''
                    if top_hit is None:
                        assign_mtb = True
                        if os.path.exists(unannotated_fasta):
                            fasta_records = SeqIO.parse(unannotated_fasta, 'fasta')
                            for mtb_record in fasta_records:
                                mtb_id = mtb_record.description.split(' ')[0]
                                mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))
                            stdout_2 = blast(unannotated_fasta, unannotated_gene_seq)
                            top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
                            if top_hit_mtb is None:
                                assign_mtb = True
                            else:
                                assign_mtb = False
                                name_to_assign = top_hit_mtb
                    else:
                        name_to_assign = top_hit
                    if assign_mtb:
                        mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                        mtb_increment = mtb_increment + 1
                        sequence_object = Seq(unannotated_gene_seq)
                        seq_record = SeqRecord(sequence_object, id=mtb_id)
                        with open(unannotated_fasta, 'a') as mtb_fasta:
                            SeqIO.write(seq_record, mtb_fasta, 'fasta')
                        name_to_assign = mtb_id
                        logger.info('Assigned new mtb id' + '\n')
                        logger.info(name_to_assign + '\n')
                    else:
                        logger.info('Assigned existing Rv or MTB' + '\n')
                        logger.info(name_to_assign + '\n')
                    update_dictionary_ltag_assignments(unannotated_gene_isolate, unannotated_gene_locus, name_to_assign)
                os.unlink(fasta_fp_to_blast)
        else:
            if len(unassigned_l_tags) > 0:
                single_gene_cluster_complete[gene] = in_dict[gene]
            continue
    with open('multi_gene_clusters.tsv', 'w') as out:
        for line in mgc_output:
            out.write(line)
    return mtb_increment, single_gene_cluster_complete


def arguments():
    parser = argparse.ArgumentParser(description='Update feature information for a set of annotations based on '
                                                 'clustering using CDHIT/MCL')
    parser.add_argument('-c', '--clusters', help='Output from cluster.py (clustered_proteins)')
    parser.add_argument('-d', '--dir', help='Directory that contains all annotations entered into the cluster.py. '
                                            'Genbank format required')
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s')
    logger = logging.getLogger(__name__)
    args = arguments()
    global isolate_update_dictionary, isolate_sequences
    isolate_update_dictionary = {}
    isolate_sequences = create_isolate_seq_dict(directory=args.dir)
    mtb_increment = find_larges_mtb_increment(annotation_directory=args.dir)
    reference_nucleotide_fastas, reference_protein_fastas = ref_seqs(args.dir)
    reference_fasta = 'ref_cdss.fasta'
    with open(reference_fasta, 'w') as ref_cds_out:
        for s in reference_protein_fastas:
            SeqIO.write(s, ref_cds_out, 'fasta')
    mtb_genes_fp = find_unannotated_genes(reference_fasta)
    pangenome_dir = args.clusters
    multi_gene_cluster, single_gene_cluster, partial_gene_cluster, candidate_novel_gene_cluster, unique_gene_cluster = \
        parse_clustered_proteins(pangenome_dir, args.dir)
    logger.info('Number of clusters with partial genes: ' + str(len(partial_gene_cluster.keys())) + '\n')
    candidate_novel_gene_cluster_complete = candidate_novel_gene_cluster.copy()
    single_gene_cluster_complete = single_gene_cluster.copy()
    mtb_increment, \
        updated_single_gene_clusters = multigene_clusters(in_dict=multi_gene_cluster,
                                                          single_gene_cluster_complete=single_gene_cluster_complete,
                                                          annotation_dir=args.dir,
                                                          unannotated_fasta=mtb_genes_fp,
                                                          mtb_increment=mtb_increment)
    single_gene_clusters(updated_single_gene_clusters)
    mtb_increment = only_ltag_clusters(in_dict=candidate_novel_gene_cluster_complete,
                                       annotation_dir=args.dir,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment,
                                       reference_fasta=reference_fasta)
    mtb_increment = singleton_clusters(singleton_dict=single_gene_cluster_complete,
                                       annotation_dir=args.dir,
                                       reference_fasta=reference_fasta,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment)
    pickle_fp = 'isolate_gene_name.pickle'
    with open(pickle_fp, 'wb') as pickle_file:
        pickle.dump(isolate_update_dictionary, pickle_file)


if __name__ == '__main__':
    main()
