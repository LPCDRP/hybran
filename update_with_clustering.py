#!/usr/bin/env python2.7

import os
import pickle
import re
import sys
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastpCommandline


GROUPHOME = os.environ['GROUPHOME']


def get_rv_sequences_from_fasta(h37rv_fasta_fp):
    rv_sequence_dict = {}
    fp = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
    rv_temp_fasta = fp.name
    fasta_records = SeqIO.parse(h37rv_fasta_fp, 'fasta')
    for record in fasta_records:
        rv_seq = str(record.seq)
        #rv = record.id
        rv_id = record.description.split('~~~')[1]
        header = '>' + rv_id + '\n'
        #print(header)
        seq = rv_seq
        rv_sequence_dict[rv_id] = seq
        fp.write(header)
        fp.write(seq)
        fp.write('\n')
    return rv_temp_fasta


def get_sequence(isolate, locus_tag):
    genbank_file = GROUPHOME + '/data/depot/annotation/' + isolate + '/annomerge/' + isolate + '.gbk'
    for record in SeqIO.parse(genbank_file, 'genbank'):
        features = record.features
        break

    for feature in features:
        if feature.type == 'CDS' and feature.qualifiers['locus_tag'][0] == locus_tag:
            seq = feature.qualifiers['translation'][0]
            return seq


def write_fasta(list_of_ids, output_name):
    list_of_records = []
    for rep in list_of_ids:
        if type(rep) is str:
            isolate_id = rep.split(',')[0]
            locus = rep.split(',')[1]
            gene = rep.split(',')[2]
        elif type(rep) is list:
            isolate_id = rep[0]
            locus = rep[1]
            gene = rep[2]

        sequence = get_sequence(isolate_id, locus)
        sequence_object = Seq(sequence)
        seq_record = SeqRecord(sequence_object, id='|'.join([isolate_id, locus, gene]))
        list_of_records.append(seq_record)

    with open(output_name, 'w') as fasta_output:
        for s in list_of_records:
            SeqIO.write(s, fasta_output, 'fasta')


def parse_clustered_proteins(pangenome_dir):
    underscore_re = re.compile('_[0-9]$')

    # Private Function #
    ##################################################################################################################
    def gff_dict():
        gff_dictionary = {}
        annotation_dir = GROUPHOME + '/data/depot/annotation/'
        isolate_fasta_list = os.listdir(GROUPHOME + '/data/genomes/')
        isolate_list = [isolate_id.split('.')[0] for isolate_id in isolate_fasta_list if isolate_id.endswith('.fasta')]
        for isolate in os.listdir(annotation_dir):
            gff = isolate + '/annomerge/' + isolate + '.gff'
            if isolate not in isolate_list and 'H37Rv' not in isolate:
                continue
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

    gffs = gff_dict()
#    pangenome_dir = GROUPHOME + '/data/depot/annotation/pangenome_04-16-18/'
    #pangenome_dir = '/home/dgunasek/projects/annomerge/pangenome_test/'
    clustered_proteins = pangenome_dir + 'clustered_proteins'
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
                #print isolate
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
#    with open('gene-clusters-dictionary.pickle', 'w') as gene_pickle:
#        pickle.dump(gene_cluster, gene_pickle)
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
            #print(line)
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
    if isolate_id not in isolate_update_dictionary.keys():
        isolate_update_dictionary[isolate_id] = {}
        isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
    else:
        isolate_dict_added = isolate_update_dictionary[isolate_id]
        if isolate_ltag in isolate_dict_added.keys():
            update_log.write('ERROR: This locus tag already has an Rv' + '\n')
            update_log.write('Locus Tag:' + '\n')
            update_log.write(isolate_id + '\n')
            update_log.write(isolate_ltag + '\n')
            update_log.write(isolate_update_dictionary[isolate_id][isolate_ltag] + '\n')
            update_log.write('Above locus tag is also annotated as: ' + new_gene_name + '\n')
        else:
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
        rep_sequence = get_sequence(rep_isolate_id, rep_locus)
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
                gene_sequence = get_sequence(gene_isolate_id, gene_locus)
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


def main():
#    print('Im HERE')
    h37rv_protein_fasta = GROUPHOME + '/resources/H37Rv-CDS-updated.fasta'
#    mtb_genes_fp = GROUPHOME + '/data/depot/hypotheticome-clinical/mtb_genes.fasta'
    mtb_genes_fp = '/home/dgunasek/projects/annomerge/mtb_genes.fasta'
    isolate_fasta_list = os.listdir(GROUPHOME + '/data/genomes/')
    global isolate_list
    isolate_list = [isolate_id.split('.')[0] for isolate_id in isolate_fasta_list if isolate_id.endswith('.fasta')]
    rv_temp_fasta = get_rv_sequences_from_fasta(h37rv_protein_fasta)
#    print(rv_temp_fasta)
    global update_log
    update_log = open('/home/dgunasek/projects/annomerge/update_with_clustering.log', 'w')
    mtb_increment = 0
    global isolate_update_dictionary
    isolate_update_dictionary = {}
    assignment_keys = {}
    pangenome_dir = GROUPHOME + '/data/depot/annotation/pangenome/'
    multi_gene_cluster, single_gene_cluster, partial_gene_cluster, candidate_novel_gene_cluster, unique_gene_cluster = \
        parse_clustered_proteins(pangenome_dir)
    #print(multi_gene_cluster.keys()[:5])
    #print(multi_gene_cluster.values()[:5])
    update_log.write('Number of clusters with partial genes: ' + str(len(partial_gene_cluster.keys())) + '\n')
    candidate_novel_gene_cluster_complete = candidate_novel_gene_cluster.copy()
    single_gene_cluster_complete = single_gene_cluster.copy()
    # 4. If a cluster has multiple H37Rv genes and L_tags, BLAST the L_tags to all genes in the cluster that are H37Rv
    # and annotate L_tag with the top hit.
    mgc_output = open('multi_gene_clusters.tsv', 'w')
    #print('Possible num multiclusters: ' + str(len(multi_gene_cluster.keys())))
    #print('Num multi genes = ' + str(len(multi_gene_cluster.keys())))
    
    num_multi = 0
    for gene in multi_gene_cluster.keys():
        num_multi += 1
        #print(num_multi) 
        #print(gene)
        output_line = ''
        unassigned_l_tags = []
        true_multi_cluster = False
        gene_elements = gene.split(',')
        for gene_in_cluster in multi_gene_cluster[gene]:
            if gene_in_cluster[1] in gene_elements or gene_in_cluster[2] in gene_elements:
                continue
            elif gene_in_cluster[1].startswith('L') and gene_in_cluster[2].startswith('L'):
                unassigned_l_tags.append(gene_in_cluster)
            else:
                true_multi_cluster = True
        #print(true_multi_cluster)
        if true_multi_cluster:
            output_line = output_line + gene
            for gene_in_cluster in multi_gene_cluster[gene]:
                gene_str = ','.join(gene_in_cluster)
                output_line = output_line + '\t' + gene_str
            output_line = output_line + '\n'
            mgc_output.write(output_line)
            all_genes_annotated = cluster_annotation_presence(multi_gene_cluster[gene])
            # If all genes in the cluster are annotated, do nothing. If there is an unannotated gene in the cluster,
            # blast it to all annotated genes in cluster and annotated with top hit
            if not all_genes_annotated:
                fasta_fp_to_blast, genes_to_annotate = get_cluster_fasta(gene, multi_gene_cluster[gene])
                for unannotated_gene in genes_to_annotate:
                    unannotated_gene_isolate = unannotated_gene[0]
                    unannotated_gene_locus = unannotated_gene[1]
                    unannotated_gene_seq = get_sequence(unannotated_gene_isolate, unannotated_gene_locus)
                    blast_command = NcbiblastpCommandline(subject=fasta_fp_to_blast,
                                                          outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
                    stdout, stderr = blast_command(stdin=unannotated_gene_seq)
                    top_hit, all_hits = identify_top_hits(stdout)
                    assign_mtb = False
                    name_to_assign = ''
                    if top_hit is None:
                        assign_mtb = True
                        if os.path.exists(mtb_genes_fp):
                            fasta_records = SeqIO.parse(mtb_genes_fp, 'fasta')
                            for mtb_record in fasta_records:
                                mtb_id = mtb_record.description.split(' ')[0]
                                mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))

                            blast_command_2 = NcbiblastpCommandline(subject=mtb_genes_fp,
                                                                    outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
                            stdout_2, stderr = blast_command_2(stdin=unannotated_gene_seq)
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
                        with open(mtb_genes_fp, 'a') as mtb_fasta:
                            SeqIO.write(seq_record, mtb_fasta, 'fasta')
                        name_to_assign = mtb_id
                        update_log.write('Assigned new mtb id' + '\n')
                        update_log.write(name_to_assign + '\n')
                    else:
                        update_log.write('Assigned existing Rv or MTB' + '\n')
                        update_log.write(name_to_assign + '\n')
                    update_dictionary_ltag_assignments(unannotated_gene_isolate, unannotated_gene_locus, name_to_assign)
                os.unlink(fasta_fp_to_blast)
        else:
            if len(unassigned_l_tags) > 0:
                #print(unassigned_l_tags)
                single_gene_cluster_complete[gene] = multi_gene_cluster[gene]
                #print(len(candidate_novel_gene_cluster_complete.keys()))
            continue
            ##### TO DO: Assign consistent gene names and Locus tags to all annotations in cluster #####
    mgc_output.close()

    # 1. Handling clusters with only one Rv/gene name in cluster
    # 1. If cluster has candidate novel genes (L_***** or L2_*****) clustered together with a H37Rv annotation, update
    #  all candidate novel genes in this cluster with the H37Rv gene name.
    rep_ltag_keys = []
    update_log.write('Number of clusters with single genes: ' + str(len(single_gene_cluster_complete.keys())) + '\n')
    for key_sgc in single_gene_cluster_complete:
        rep_ltag = False
        key_elements_sgc = key_sgc.split(',')
        if key_elements_sgc[1].startswith('Rv'):
            gene_to_add = key_elements_sgc[1]
        elif key_elements_sgc[2].startswith('Rv'):
            gene_to_add = key_elements_sgc[2]
#        elif (key_elements_sgc[1].startswith('L_') and key_elements_sgc[2].startswith('L_')) or \
#                (key_elements_sgc[1].startswith('L2_') and key_elements_sgc[2].startswith('L2_')):
        elif (key_elements_sgc[1].startswith('L') and key_elements_sgc[2].startswith('L')):
            rep_ltag = True
            rep_ltag_keys.append(key_sgc)
        else:
            gene_to_add = key_elements_sgc[2]
        # If representative sequence is an Rv and not an L-tag
        if not rep_ltag:
            genes_in_cluster = single_gene_cluster_complete[key_sgc]
            for gene in genes_in_cluster:
#                if (gene[1].startswith('L_') and gene[2].startswith('L_')) or \
#                        (gene[1].startswith('L2_') and gene[2].startswith('L2_')):
                if (gene[1].startswith('L') and gene[2].startswith('L')):
                    update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                else:
                    continue
        # If representative sequence is a L-tag
        else:
            genes_in_cluster = single_gene_cluster_complete[key_sgc]
            gene_to_add = ''
            for gene in genes_in_cluster:
#                if (gene[1].startswith('L_') and gene[2].startswith('L_')) or \
#                        (gene[1].startswith('L2_') and gene[2].startswith('L2_')):
                if (gene[1].startswith('L') and gene[2].startswith('L')):
                    continue
                else:
                    if not gene[1].startswith('L'):
                        gene_to_add = gene[1]
                    else:
                        gene_to_add = gene[2]
                    break
            if len(gene_to_add) == 0:
                update_log.write('No Rv found in this cluster' + '\n')
                update_log.write(key_sgc + '\n')
                update_log.write(genes_in_cluster + '\n')
            else:
                for gene in genes_in_cluster:
#                    if (gene[1].startswith('L_') and gene[2].startswith('L_')) or \
#                            (gene[1].startswith('L2_') and gene[2].startswith('L2_')):
                    if (gene[1].startswith('L') and gene[2].startswith('L')):
                        update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                    else:
                        continue


    # 2. If a cluster has only candidate novel genes (with no gene names), then BLAST representative sequence to H37Rv
    # and if there is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster
    # is annotated with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes
    # in the cluster.
    update_log.write('Number of clusters with only L-tags genes: ' + str(len(candidate_novel_gene_cluster_complete.keys())) +
                     '\n')
    for rep_gene in candidate_novel_gene_cluster_complete.keys():
        l_tag_verify_fp = '/home/dgunasek/projects/annomerge/verify_ltag_clusters.tsv'
        rep_isolate_id = rep_gene.split(',')[0]
        rep_locus = rep_gene.split(',')[1]
        rep_gene_name = rep_gene.split(',')[2]
        rep_sequence = get_sequence(rep_isolate_id, rep_locus)
        # Writing representative amino acid sequence to a temp file to check if all genes in cluster share identity with
        #  this sequence
        rep_fp = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
        rep_temp_fasta = rep_fp.name
        rep_sequence_object = Seq(rep_sequence)
        rep_record = SeqRecord(rep_sequence_object, id='|'.join([rep_isolate_id, rep_locus, rep_gene_name]))
        with open(rep_temp_fasta, 'w') as fasta_tmp:
            SeqIO.write(rep_record, fasta_tmp, 'fasta') ###
        # Blasting representative sequence with H37Rv
        blast_command = NcbiblastpCommandline(subject=rv_temp_fasta,
                                              outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
        stdout, stderr = blast_command(stdin=rep_sequence)
        top_hit, all_hits = identify_top_hits(stdout)
        assign_mtb = False
        name_to_assign = ''
        if top_hit is None:
            assign_mtb = True
            if os.path.exists(mtb_genes_fp):
                fasta_records = SeqIO.parse(mtb_genes_fp, 'fasta')
                for mtb_record in fasta_records:
                    mtb_id = mtb_record.description.split(' ')[0]
                    mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))

                blast_command_2 = NcbiblastpCommandline(subject=mtb_genes_fp,
                                                        outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
                stdout_2, stderr = blast_command_2(stdin=rep_sequence)
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
            with open(mtb_genes_fp, 'a') as mtb_fasta:
                SeqIO.write(seq_record, mtb_fasta, 'fasta')
            name_to_assign = mtb_id
            update_log.write('Assigned new mtb id' + '\n')
            update_log.write(name_to_assign + '\n')
        else:
            update_log.write('Assigned existing Rv or MTB' + '\n')
            update_log.write(name_to_assign + '\n')
        for gene_in_cluster in candidate_novel_gene_cluster_complete[rep_gene]:
#            gene_in_cluster_seq = get_sequence(gene_in_cluster[0], gene_in_cluster[1])
#            blast_to_rep = NcbiblastpCommandline(subject=rep_temp_fasta,
#                                                 outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
#            stdout_rep, stderr_rep = blast_to_rep(stdin=gene_in_cluster_seq)
#            top_hit_rep, all_hits_rep = identify_top_hits(stdout_rep)
#            if top_hit_rep is None:
#                update_log.write('Gene in cluster does not pass threshold of identity with representative sequence' +
#                                 '\n')
#                update_log.write(rep_isolate_id + '\n')
#                update_log.write(rep_locus + '\n')
#                update_log.write(rep_sequence + '\n')
#                update_log.write(rep_temp_fasta + '\n')
#                update_log.write(gene_in_cluster[0] + '\n')
#                update_log.write(gene_in_cluster[1] + '\n')
#                update_log.write(gene_in_cluster_seq + '\n')
#                rep_seq_id = '|'.join([rep_isolate_id, rep_locus])
#                gene_seq_id = '|'.join([gene_in_cluster[0], gene_in_cluster[1]])
#                with open(l_tag_verify_fp, 'a+') as verify_ltag:
#                    verify_ltag.write(rep_seq_id + '\t' + gene_seq_id + '\n')
            update_dictionary_ltag_assignments(gene_in_cluster[0], gene_in_cluster[1], name_to_assign)
            #else:
            #    update_dictionary_ltag_assignments(gene_in_cluster[0], gene_in_cluster[1], name_to_assign)
        os.unlink(rep_temp_fasta)


    # 3. If a cluster has only one gene (with no gene names), then BLAST representative sequence to H37Rv and if there
    # is a hit with specified amino acid and coverage thresholds, all candidate novel genes in the cluster is annotated
    #  with the H37Rv gene. If the representative does not hit a H37Rv, assign a MTB locus tag to the genes in the
    # cluster.
    update_log.write('Number of singleton clusters with single genes: ' + str(len(unique_gene_cluster)) + '\n')
    for single_gene in unique_gene_cluster:
        isolate_id = single_gene[0]
        locus_tag = single_gene[1]
        gene_name = single_gene[2]
        if gene_name.startswith('L') and locus_tag.startswith('L'):
            gene_sequence = get_sequence(isolate_id, locus_tag)
            blast_command = NcbiblastpCommandline(subject=rv_temp_fasta,
                                                  outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
            stdout, stderr = blast_command(stdin=gene_sequence)
            top_hit, all_hits = identify_top_hits(stdout)
            assign_mtb = False
            name_to_assign = ''
            if top_hit is None:
                assign_mtb = True
                if os.path.exists(mtb_genes_fp):
                    fasta_records = SeqIO.parse(mtb_genes_fp, 'fasta')
                    for mtb_record in fasta_records:
                        mtb_id = mtb_record.description.split(' ')[0]
                        mtb_increment = max(mtb_increment, int(mtb_id.split('MTB')[1]))

                    blast_command_2 = NcbiblastpCommandline(subject=mtb_genes_fp,
                                                            outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
                    stdout_2, stderr = blast_command_2(stdin=gene_sequence)
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
                with open(mtb_genes_fp, 'a') as mtb_fasta:
                    SeqIO.write(seq_record, mtb_fasta, 'fasta')
                name_to_assign = mtb_id
                update_log.write('Assigned new mtb id' + '\n')
                update_log.write(name_to_assign + '\n')
            else:
                update_log.write('Assigned existing Rv or MTB' + '\n')
                update_log.write(name_to_assign + '\n')
            update_dictionary_ltag_assignments(isolate_id, locus_tag, name_to_assign)
        else:
            continue

    update_log.close()
    os.unlink(rv_temp_fasta)
    pickle_fp = '/home/dgunasek/projects/annomerge/isolate_gene_name.pickle'
    with open(pickle_fp, 'wb') as pickle_file:
        pickle.dump(isolate_update_dictionary, pickle_file)


if __name__ == '__main__':
    main()
