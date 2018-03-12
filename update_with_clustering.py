#!/usr/bin/env python2.7

import os
import pickle
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

GROUPHOME = os.environ['GROUPHOME']
underscore_re = re.compile('_[0-9]$')

def parse_clustered_proteins():

    # Private Function #
    ##################################################################################################################
    def gff_dict():
        gff_dictionary = {}
        annotation_dir = GROUPHOME + '/data/depot/annotation/'
        for isolate in os.listdir(annotation_dir):
            gff = isolate + '/annomerge/' + isolate + '.gff'
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

    def write_fasta(list_of_ids):
        def get_sequence(isolate, locus_tag):
            genbank_file = GROUPHOME + '/data/depot/annotation/' + isolate + '/annomerge/' + isolate + '.gbk'
            features = SeqIO.read(genbank_file, 'genbank').features
            for feature in features:
                if feature.type == 'CDS' and features.qualifiers['locus_tag'][0] == locus_tag:
                    seq = feature.qualifiers['translation']
                    return seq
        list_of_records = []
        for rep in list_of_ids:
            isolate_id = rep.split(',')[0]
            locus = rep.split(',')[1]
            gene = rep.split(',')[2]
            sequence = get_sequence(isolate_id, locus)
            seq_record = SeqRecord(sequence, id=isolate_id, description='|'.join([locus, gene]))
            list_of_records.append(seq_record)
        with open('representative-clusters.fasta', 'w') as fasta_output:
            for s in list_of_records:
                SeqIO.write(s, fasta_output, 'fasta')


    ##################################################################################################################

    gffs = gff_dict()
    pangenome_dir = GROUPHOME + '/data/depot/annotation/pangenome/'
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


if __name__ == '__main__':
    parse_clustered_proteins()
