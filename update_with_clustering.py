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
from sets import Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastpCommandline


def parse_clustered_proteins(clustered_proteins, annotations):
    underscore_re = re.compile('_[0-9]$')

    # Private Function #
    ##################################################################################################################
    def gff_dict(annotation_dir):
        gff_dictionary = {}
        for gff in os.listdir(annotation_dir):
            isolate_id_ltag = {}
            try:
                isolate_id = gff.split('.')[0]
                if gff.endswith('.gff'):
                    with open(annotation_dir + gff, 'r') as gff_file:
                        for line in gff_file:
                            if line.startswith('#'):
                                continue
                            gene = ''
                            column = line.rstrip('\n').split('\t')
                            if len(column) >= 8 and 'CDS' in column[2]:
                                info = column[8].split(';')
                                gene_id = ''.join([i.split('=')[1] for i in info if i.startswith('ID=')])
                                locus_tag = ''.join([i.split('=')[1] for i in info if i.startswith('locus_tag=')])
                                gene = ','.join([i.split('=')[1] for i in info if i.startswith('gene=')])
                                if not locus_tag.startswith('Rv') and not locus_tag.startswith('L'):
                                    gene = locus_tag
                                    locus_tag = ''
                                isolate_id_ltag[gene_id] = (locus_tag, gene)
                        gff_dictionary[isolate_id] = isolate_id_ltag
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
            if (len(isolates_ids) == 1 and
                    (representative_ltag_gene_tup[1].startswith('L')) or
                     underscore_re.search(representative_ltag_gene_tup[1])):
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
                    any([gene[1].startswith('L') for gene in cluster_list]):
                different_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            #
            elif (len(unique_genes) == 1 or len(unique_locus) == 1) and \
                    any([gene[1].startswith('L') for gene in cluster_list]):
                same_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            # L tag only clusters
            elif all(locus[1].startswith('L') for locus in cluster_list):
                l_tag_only_clusters[representative] = cluster_list_w_isolate
    return [different_genes_cluster_w_ltags, same_genes_cluster_w_ltags, l_tag_only_clusters,
            unique_genes_list]


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
        if line_elements and line_elements != [' ']:
            iden = float(line_elements[5])
            qcov = (float(line_elements[4]) / float(line_elements[1])) * 100.0
            scov = (float(line_elements[4]) / float(line_elements[3])) * 100.0
            if iden >= identity and qcov >= coverage and scov >= coverage and not line_elements[2].startswith('L'):
                rv_hit = True
                if '|' in line_elements[2]:
                    corresponding_rv_hit = line_elements[2].split('|')[-1]
                else:
                    corresponding_rv_hit = line_elements[2]
                all_hits_dict[corresponding_rv_hit] = iden
            else:
                continue
    if rv_hit:
        top_hit = get_top_hit(all_hits_dict)
    else:
        top_hit = None
    return top_hit, all_hits_dict


def update_dictionary_ltag_assignments(isolate_id, isolate_ltag, new_gene_name):
    logger = logging.getLogger('NameGene')
    if isolate_id not in isolate_update_dictionary.keys():
        isolate_update_dictionary[isolate_id] = {}
        isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
        logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
    else:
        isolate_dict_added = isolate_update_dictionary[isolate_id]
        if isolate_ltag not in isolate_dict_added.keys():
            logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
            isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
    return


def get_cluster_fasta(rep, cluster_list):
    rep_isolate_id = rep.split(',')[0]
    rep_locus = rep.split(',')[1]
    rep_gene_name = rep.split(',')[2]
    cluster_fasta_fp = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
    cluster_temp_fasta = cluster_fasta_fp.name
    added_seq = []
    unannotated_genes_list = []
    fasta_tmp = open(cluster_temp_fasta, 'w')
    if not (rep_locus.startswith('L') and rep_gene_name.startswith('L')):
        rep_sequence = isolate_sequences[rep_isolate_id][rep_locus]
        rep_sequence_object = rep_sequence
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
                gene_sequence = isolate_sequences[gene_isolate_id][gene_locus]
                gene_sequence_object = gene_sequence
                added_seq.append(gene_id)
                gene_record = SeqRecord(gene_sequence_object, id=gene_id)
                SeqIO.write(gene_record, fasta_tmp, 'fasta')
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


def find_largest_mtb_increment(unannotated_fasta):
    mtbs = []
    for record in SeqIO.parse(unannotated_fasta, 'fasta'):
        mtbs.append(record.id)
    largest_mtb = sorted(mtbs)[-1]
    if largest_mtb:
        return int(largest_mtb.replace('MTB', ''))
    return 0


def ref_seqs(gbk_dir):

    def create_allseq_dict(fa):
        """Parses FASTA input, creates dictionary with sequence ID as key and the sequence as the value"""
        fasta_dict = {}
        for record in SeqIO.parse(fa, "fasta"):
            sequence = str(record.seq)
            fasta_dict[record.id] = sequence
        return fasta_dict

    def run_cdhit(nproc, input, output):
        """Runs cdhit with a sequence identity threshold of 0.98."""
        cmd = ['cdhit',
               '-i', input,
               '-o', output,
               '-c', '0.99',
               '-T', str(nproc),
               '-g', '1',
               '-s', '1.0',
               '-d', '256',
               '-A', '1.0']
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        return out

    def create_reps_dict(in_clusters):
        """Creates a dictionary with the representative sequence ID from each cluster as the key
        and the corresponding sequences as the value."""
        clusters = []
        rep_dict = {}
        rep = ''
        reps = []
        with open(in_clusters, 'r') as f:
            for line in f:
                # At new cluster
                if line.startswith('>'):
                    if not rep:
                        continue
                    rep_dict[rep] = clusters
                    clusters = []
                    continue
                # Representatives of cluster
                else:
                    if line.rstrip('\n').endswith('*'):
                        rep = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                        reps.append(rep)
                    else:
                        gene = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                        clusters.append(gene)
        return rep_dict, reps

    def create_reps_fasta(output, reps, rep_seq_id_dict):
        """Creates a FASTA file that has sequence ID and sequence for each representative"""
        seq_list = []
        for id_key in reps:
            record = SeqRecord(Seq(rep_seq_id_dict[id_key]), id=id_key, description="")
            seq_list.append(record)
        SeqIO.write(seq_list, output, "fasta")

    def cd_hit(nproc, fasta, out):
        logger = logging.getLogger('CDHIT')
        logger.debug('Running CDHIT on reference annotations')
        cdhit_stdout = run_cdhit(nproc, fasta, out)
        OGdict = create_allseq_dict(fasta)
        REPdict, rep_list = create_reps_dict(out + ".clstr")
        create_reps_fasta(out, rep_list, OGdict)

    def grep_seqs(gff):
        cmd = ['grep', 'translation=', gff]
        translations = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        return [line for line in translations.stdout]
    protein_cds = []
    isolate_seqs = {}
    for gff in os.listdir(gbk_dir):
        if gff.endswith('.gff'):
            raw_out = grep_seqs(gbk_dir + gff)
            gff_name = gff.split('.')[0]
            isolate_seqs[gff_name] = {}
            for line in raw_out:
                if line.split('\t')[2] == 'CDS':
                    locus_tag = [i.split('=')[1].rstrip('\n') for i in line.split(';') if i.startswith('locus_tag=')][0]
                    gene = [i.split('=')[1].rstrip('\n') for i in line.split(';') if i.startswith('gene=')][0]
                    translation = [i.split('=')[1] for i in line.split(';') if i.startswith('translation=')][0]
                    record = SeqRecord(Seq(translation.rstrip('\n')),
                                       id=gene,
                                       description=gff_name)
                    isolate_seqs[gff_name][locus_tag] = Seq(translation.rstrip('\n'))
                    protein_cds.append(record)
    all_proteins = 'ref_cdss_protein-all.fasta'
    final_proteins = 'ref_cdss_protein.fasta'
    with open(all_proteins, 'w') as ref_cds_out:
        for s in protein_cds:
            SeqIO.write(s, ref_cds_out, 'fasta')
    cd_hit(1, all_proteins, final_proteins)
    return final_proteins, isolate_seqs


def blast(subject, stdin_seq):
    blast_command = NcbiblastpCommandline(subject=subject,
                                          outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
    stdout, stderr = blast_command(stdin=str(stdin_seq))
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


def singleton_clusters(singleton_dict, reference_fasta, unannotated_fasta, mtb_increment):
    logger = logging.getLogger('SingletonClusters')
    # 3. If a cluster has only one gene (with no gene names), then BLAST representative sequence to H37Rv and if there
    # is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster is annotated
    #  with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes in the
    # cluster.
    logger.debug('Number of singleton clusters with single genes: ' + str(len(singleton_dict)))
    for single_gene in singleton_dict:
        isolate_id = single_gene[0]
        locus_tag = single_gene[1]
        gene_name = single_gene[2]
        if gene_name.startswith('L') and locus_tag.startswith('L'):
            gene_sequence = isolate_sequences[isolate_id][locus_tag]
            stdout = blast(reference_fasta, gene_sequence)
            top_hit, all_hits = identify_top_hits(stdout)
            assign_mtb = False
            name_to_assign = top_hit
            if top_hit is None:
                assign_mtb = True
                stdout_2 = blast(unannotated_fasta, gene_sequence)
                top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
                if top_hit_mtb is not None:
                    assign_mtb = False
                    name_to_assign = top_hit_mtb
            if assign_mtb:
                mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                mtb_increment = mtb_increment + 1
                seq_record = SeqRecord(gene_sequence, id=mtb_id)
                with open(unannotated_fasta, 'a') as mtb_fasta:
                    SeqIO.write(seq_record, mtb_fasta, 'fasta')
                name_to_assign = mtb_id
            update_dictionary_ltag_assignments(isolate_id, locus_tag, name_to_assign)
        else:
            continue
    return mtb_increment


def only_ltag_clusters(in_dict, reference_fasta, unannotated_fasta, mtb_increment):
    logger = logging.getLogger('NewGeneClusters')
    # 2. If a cluster has only candidate novel genes (with no gene names), then BLAST representative sequence to H37Rv
    # and if there is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster
    # is annotated with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes
    # in the cluster.
    logger.debug('Number of clusters with no gene names: ' + str(len(in_dict.keys())))
    new_unannotated_genes = []
    rep_records = []
    for rep_gene in in_dict.keys():
        rep_isolate_id = rep_gene.split(',')[0]
        rep_locus = rep_gene.split(',')[1]
        rep_gene_name = rep_gene.split(',')[2]
        rep_sequence = isolate_sequences[rep_isolate_id][rep_locus]
        rep_record = SeqRecord(rep_sequence, id='|'.join([rep_isolate_id, rep_locus, rep_gene_name]),
                               description='')
        rep_records.append(rep_record)

    rep_fp = tempfile.NamedTemporaryFile(suffix='.fasta', dir='/tmp/')
    rep_temp_fasta = rep_fp.name
    with open(rep_temp_fasta, 'w') as fasta_tmp:
        for s in rep_records:
            SeqIO.write(s, fasta_tmp, 'fasta')
    blastcmd = NcbiblastpCommandline(subject=reference_fasta,
                                     query=rep_temp_fasta,
                                     outfmt='"6 qseqid qlen sseqid slen qlen length pident qcovs"')
    stdout, stderr = blastcmd()
    rep_fp.close()
    hits = {}
    for line in stdout.split('\n'):
        if line:
            try:
                hits[line.split('\t')[0]] += [line]
            except KeyError:
                hits[line.split('\t')[0]] = [line]
    for qid, lines in hits.iteritems():
        top_hit, all_hits = identify_top_hits('\n'.join(lines))
        isolate = qid.split('|')[0]
        locus = qid.split('|')[1]
        assign_mtb = False
        name_to_assign = top_hit
        if top_hit is None:
            assign_mtb = True
            rep_seq = isolate_sequences[qid.split('\t')[0].split('|')[0]][qid.split('\t')[0].split('|')[1]]
            stdout_2 = blast(unannotated_fasta, rep_seq)
            top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
            if top_hit_mtb is not None:
                assign_mtb = False
                name_to_assign = top_hit_mtb
        if assign_mtb:
            mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
            mtb_increment = mtb_increment + 1
            rep_sequence = isolate_sequences[qid.split('\t')[0].split('|')[0]][qid.split('\t')[0].split('|')[1]]
            seq_record = SeqRecord(rep_sequence, id=mtb_id)
            new_unannotated_genes.append(seq_record)
            name_to_assign = mtb_id
        update_dictionary_ltag_assignments(isolate, locus, name_to_assign)
    with open(unannotated_fasta, 'a') as mtb_fasta:
        for s in new_unannotated_genes:
            SeqIO.write(s, mtb_fasta, 'fasta')
    return mtb_increment


def single_gene_clusters(single_gene_dict):
    logger = logging.getLogger('ClustersWithOneGeneName')
    # 1. Handling clusters with only one Rv/gene name in cluster
    # 1. If cluster has candidate novel genes (L_***** or L2_*****) clustered together with a H37Rv annotation, update
    #  all candidate novel genes in this cluster with the H37Rv gene name.
    rep_ltag_keys = []
    logger.debug('Number of clusters with single genes: ' + str(len(single_gene_dict.keys())))
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
            if len(gene_to_add) != 0:
                for gene in genes_in_cluster:
                    if gene[1].startswith('L') and gene[2].startswith('L'):
                        update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                    else:
                        continue


def multigene_clusters(in_dict, single_gene_cluster_complete, unannotated_fasta, mtb_increment):
    # 4. If a cluster has multiple H37Rv genes and L_tags, BLAST the L_tags to all genes in the cluster that are H37Rv
    # and annotate L_tag with the top hit.
    logger = logging.getLogger('ClustersWithManyGeneNames')
    logger.debug('Number of clusters that have many gene names and genes with no name ' + str(len(in_dict.keys())))
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
                fasta_fp_to_blast, genes_to_annotate = get_cluster_fasta(gene, in_dict[gene])
                for unannotated_gene in genes_to_annotate:
                    unannotated_gene_isolate = unannotated_gene[0]
                    unannotated_gene_locus = unannotated_gene[1]
                    unannotated_gene_seq = isolate_sequences[unannotated_gene_isolate][unannotated_gene_locus]
                    stdout = blast(fasta_fp_to_blast, unannotated_gene_seq)
                    top_hit, all_hits = identify_top_hits(stdout)
                    assign_mtb = False
                    name_to_assign = top_hit
                    if top_hit is None:
                        assign_mtb = True
                        stdout_2 = blast(unannotated_fasta, unannotated_gene_seq)
                        top_hit_mtb, all_hits_mtb = identify_top_hits(stdout_2)
                        if top_hit_mtb is not None:
                            assign_mtb = False
                            name_to_assign = top_hit_mtb
                    if assign_mtb:
                        mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                        mtb_increment = mtb_increment + 1
                        seq_record = SeqRecord(unannotated_gene_seq, id=mtb_id)
                        with open(unannotated_fasta, 'a') as mtb_fasta:
                            SeqIO.write(seq_record, mtb_fasta, 'fasta')
                        name_to_assign = mtb_id
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


def add_gene_names_to_gbk(mtb_pickle, gbk_dir, suffix):
    logger = logging.getLogger('UpdateGenbank')
    gbk_files = os.listdir(gbk_dir)
    isolates = [isolate_id.split('.')[0] for isolate_id in gbk_files if isolate_id.endswith('.gbk')]
    for isolate in mtb_pickle.keys():
        if isolate not in isolates:
            logger.error('Isolate ' + isolate + ' absent in ' + gbk_dir)
        else:
            genbank_file = gbk_dir + '/' + isolate + '.gbk'
            isolate_records = list(SeqIO.parse(genbank_file, 'genbank'))
            update_mtb_dict = mtb_pickle[isolate]
            locus_to_update = update_mtb_dict.keys()
            modified_locus = []
            rec_num = 1
            for rec in isolate_records:
                if rec_num > 1:
                    break
                rec_num += 1
                for feature in rec.features:
                    if feature.type != 'CDS':
                        continue
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    gene_name = feature.qualifiers['gene'][0]
                    try:
                        gene_synonym = feature.qualifiers['gene_synonym']
                    except KeyError:
                        gene_synonym = []
                    if locus_tag in locus_to_update and \
                            locus_tag not in modified_locus:
                        if gene_name.startswith('L'):
                            feature.qualifiers['gene'] = update_mtb_dict[locus_tag]
                        elif gene_name == update_mtb_dict[locus_tag]:
                            modified_locus.append(locus_tag)
                            continue
                        else:
                            logger.debug('Discordant assignment of gene name')
                            logger.debug('Original gene name: ' + gene_name)
                            logger.debug('New gene name: ' + update_mtb_dict[locus_tag])
                            if 'gene_synonym' in feature.qualifiers.keys():
                                gene_synonym.append(update_mtb_dict[locus_tag])
                            else:
                                feature.qualifiers['gene_synonym'] = [update_mtb_dict[locus_tag]]
                        modified_locus.append(locus_tag)
                    elif locus_tag in locus_to_update and locus_tag in modified_locus:
                        logger.debug('Updated locus tag previously')
                    else:
                        continue
                if len(Set(locus_to_update).intersection(Set(modified_locus))) < len(locus_to_update):
                    logger.warning('The following locus_tags are missing in the genbank file')
                    logger.warning(Set(locus_to_update).difference(Set(modified_locus)))
            if not suffix:
                SeqIO.write(isolate_records, genbank_file, 'genbank')
            else:
                SeqIO.write(isolate_records, gbk_dir + '/' + isolate + suffix + '.gbk', 'genbank')


def arguments():
    parser = argparse.ArgumentParser(description='Update feature information for a set of annotations based on '
                                                 'clustering using CDHIT/MCL. If the Genbank files should not be '
                                                 'overwritten, pass in -s/--suffix [suffix-name] to create new '
                                                 'Genbank files.')
    parser.add_argument('-c', '--clusters', help='Output from cluster.py (clustered_proteins)')
    parser.add_argument('-d', '--dir', help='Directory that contains all annotations entered into the cluster.py. '
                                            'Genbank format required')
    parser.add_argument('-s', '--suffix', help='Suffix output file name if the input Genbank files should not be '
                                               'overwritten. Default is to overwrite the existing Genbank files',
                        required=False)
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    return parser.parse_args()


def main():
    args = arguments()
    logger = logging.getLogger('UpdateAnnotations')
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s:INFO:%(name)s:%(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:INFO:%(name)s:%(message)s')

    global isolate_update_dictionary, isolate_sequences
    isolate_update_dictionary = {}
    logger.info('Retrieving reference protein sequences from Genbank files in ' + args.dir)
    reference_protein_fastas, isolate_sequences = ref_seqs(gbk_dir=args.dir)
    logger.info('Identifying unannotated proteins (signified by absence of a gene name)')
    mtb_genes_fp = find_unannotated_genes(reference_protein_fasta=reference_protein_fastas)
    mtb_increment = find_largest_mtb_increment(unannotated_fasta=mtb_genes_fp)
    logger.info('Parsing ' + args.clusters)
    clusters = parse_clustered_proteins(clustered_proteins=args.clusters,
                                        annotations=args.dir)
    multi_gene_cluster = clusters[0]
    single_gene_cluster = clusters[1]
    candidate_novel_gene_cluster = clusters[2]
    unique_gene_cluster = clusters[3]
    candidate_novel_gene_cluster_complete = candidate_novel_gene_cluster.copy()
    single_gene_cluster_complete = single_gene_cluster.copy()
    logger.info('Annotating genes that do not have a gene name but cluster with genes that have many gene names')
    mtb_increment, \
        updated_single_gene_clusters = multigene_clusters(in_dict=multi_gene_cluster,
                                                          single_gene_cluster_complete=single_gene_cluster_complete,
                                                          unannotated_fasta=mtb_genes_fp,
                                                          mtb_increment=mtb_increment)
    logger.info('Annotating genes that do not have a gene name but cluster with genes that have the same gene name')
    single_gene_clusters(updated_single_gene_clusters)
    mtb_increment = only_ltag_clusters(in_dict=candidate_novel_gene_cluster_complete,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment,
                                       reference_fasta=reference_protein_fastas)
    logger.info('Annotating genes that do not have a gene name and exist in a single isolate')
    mtb_increment = singleton_clusters(singleton_dict=unique_gene_cluster,
                                       reference_fasta=reference_protein_fastas,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment)
    logger.info('Updating Genbank files in ' + args.dir)
    add_gene_names_to_gbk(mtb_pickle=isolate_update_dictionary,
                          gbk_dir=args.dir,
                          suffix=args.suffix)


if __name__ == '__main__':
    main()
