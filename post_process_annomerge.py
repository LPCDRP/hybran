#!/usr/bin/env python2.7


__author__ = "Deepika Gunasekaran"
__author__ = "Sarah Ramirez-Busby"
__version__ = "0.0.1"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@sdsu.edu"
__status__ = "Development"


import tempfile
from Bio import SeqIO
import os
import collections
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline


def get_rv_gene_positions():
    genes_positions = {}
    genes_tsv_lines = open(genes_tsv_fp, 'r').readlines()
    for line in genes_tsv_lines[1:]:
        line_elements = line.strip().split('\t')
        if line_elements[2] != 'CDS':
            continue
        genes_positions[(line_elements[0])] = (int(line_elements[8]), int(line_elements[9]))
    return genes_positions


GROUPHOME = os.environ['GROUPHOME']
input_dir = GROUPHOME + '/data/annotation/backup-2018-11-13/'
output_dir = GROUPHOME + '/data/pangenome/gff-copies/all/temp-gbks-processed/'
all_files = os.listdir(input_dir)
gbk_files = [i for i in all_files if i.endswith('.gbk')]
h37rv_fasta_fp = GROUPHOME + '/data/genomes/H37Rv-NCBI.fasta'
h37rv_records = list(SeqIO.parse(h37rv_fasta_fp, "fasta"))
h37rv_sequence = h37rv_records[0].seq
h37rv = h37rv_sequence
genes_tsv_fp = GROUPHOME + '/resources/mtb-reconstruction/genes.tsv'
rv_genes_positions = get_rv_gene_positions()


def get_isolate_gene_seq(isolate, coordinates):
    ### AUTHOR: Sarah Ramirez-Busby ###
    ### MODIFIED BY: Deepika Gunasekaran ###
    for record in SeqIO.parse(GROUPHOME + '/data/genomes/' + isolate + '.fasta', 'fasta'):
        nuc = record.seq
        break
    return nuc[coordinates[0]:coordinates[1]]


def pick_best_hit(isolate, ratt_feature, prokka_feature):
    ### AUTHOR: Sarah Ramirez-Busby ###
    ### MODIFIED BY: Deepika Gunasekaran ###
    gene = ratt_feature.qualifiers['locus_tag'][0]
    if gene not in rv_genes_positions.keys():
        print('This is a merged gene. Therefore, RATT annotation is kept')
        print(ratt_feature)
        take_ratt = True
        return take_ratt
    #startpos, stoppos = get_rv_seq(gene)
    startpos = rv_genes_positions[gene][0]
    stoppos = rv_genes_positions[gene][1]
    prokka_seq = get_isolate_gene_seq(isolate, (prokka_feature.location.start, prokka_feature.location.end))
    ratt_seq = get_isolate_gene_seq(isolate, (ratt_feature.location.start, ratt_feature.location.end))
    h37rv_gene_seq = h37rv[startpos - 1: stoppos]
    with open(gene + '.fasta', 'w') as f:
        SeqIO.write(SeqRecord(h37rv_gene_seq, id=gene, description=''), f, 'fasta')
    #blast_out = blast(gene + '.fasta', prokka_seq)
    #print(ratt_seq)
    blast_out = blast(gene + '.fasta', ratt_seq)
    os.remove(gene + '.fasta')
    blast_results = False
    for i in blast_out:
        if i.startswith('Query'):
            sstart = int(i.split('\t')[8])
            send = int(i.split('\t')[9])
            blast_results = True
    if blast_results:
        # If the gene is positive strand
        if ratt_feature.location.strand == '+':
            # The start of the prokka gene minus where the prokka sequence aligns to H37Rv plus one to include the position
            #prokka_prom_end = prokka_feature.location.start.position - sstart + 1
            #prokka_prom_start = prokka_prom_end - 40
            ratt_prom_end = ratt_feature.location.start.position - sstart + 1
            ratt_prom_start = ratt_prom_end - 40
            rv_prom_end = startpos - 40
            rv_seq = h37rv[startpos - 1: rv_prom_end]
        # If the gene is negative strand
        else:
            rv_prom_end = stoppos + 40
            # Where the prokka sequence aligns to the H37Rv seq plus the end coordinate of the prokka gene
            #prokka_prom_start = send + prokka_feature.location.end.position
            #prokka_prom_end = prokka_prom_start + 40
            ratt_prom_start = send + ratt_feature.location.end.position
            ratt_prom_end = ratt_prom_start + 40
            rv_seq = h37rv[stoppos - 1: rv_prom_end]
        with open(gene + '-prom.fasta', 'w') as f:
            SeqIO.write(SeqRecord(rv_seq, id=gene, description=''), f, 'fasta')
        take_ratt = True
        #prom_blast = blast(gene + '-prom.fasta', get_isolate_gene_seq(isolate, (prokka_prom_start, prokka_prom_end)))
        prom_blast = blast(gene + '-prom.fasta', get_isolate_gene_seq(isolate, (ratt_prom_start, ratt_prom_end)))
        os.remove(gene + '-prom.fasta')
        #print prom_blast
        for i in prom_blast:
            if i.startswith('Query'):
                if int(i.split('\t')[2]) < 100.0:
                    take_ratt = False
                else:
                    take_ratt = True
    else:
        print('No Blast Hits for the RATT annotation, hence Prokka annotation is chosen: ' + gene)
        take_ratt = False
        print('RATT seq')
        print(ratt_seq)
        print('Prokka Seq')
        print(prokka_seq)
        print(blast_out) 
        print('H37Rv Seq')
        print(h37rv_gene_seq)
        print(ratt_feature)
    return take_ratt


def blast(seq1, seq2):
    ### AUTHOR: Sarah Ramirez-Busby ###
    blast_to_rv = NcbiblastnCommandline(subject=seq1,
                                        outfmt='"7 qseqid sseqid pident length mismatch '
                                               'gapopen qstart qend sstart send evalue bitscore"')
    stdout, stderr = blast_to_rv(str(seq2))
    stdout_elements = stdout.split('\n')
    #print(stdout_elements)
    return stdout_elements


def get_ordered_features(feature_list):
    """
    :param feature_list: list of features of type SeqFeature
    :return: sorted list of features, sorted by genomic location
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



def get_prom_for_rv(feature_list, source_seq):
    rv_prom_dict = {}
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
        fp_prom = tempfile.NamedTemporaryFile(suffix='_prom.fasta', delete=False, mode='w')
        rv_prom_dict[locus] = fp_prom.name
        header_prom = '>' + locus + '_prom\n'
        seq_prom = prom_seq + '\n'
        fp_prom.write(header_prom)
        fp_prom.write(seq_prom)
        fp_prom.close()
    return rv_prom_dict


def main():
    h37rv_embl_fp = GROUPHOME + '/bin/ratt-code-18/H37Rv-NC_TSS_000962.3.embl'
    h37rv_embl_record = SeqIO.read(h37rv_embl_fp, 'embl')
    h37rv_features = h37rv_embl_record.features
    h37rv_prom_fp_dict = get_prom_for_rv(h37rv_features, h37rv_sequence)
    for isolate_gbk in gbk_files:
        output_fp = output_dir + isolate_gbk
        isolate_id = isolate_gbk.split('.')[0]
        print(isolate_id)
        isolate_fasta_fp = GROUPHOME + '/data/genomes/' + isolate_id + '.fasta'
        isolate_records = list(SeqIO.parse(isolate_fasta_fp, "fasta"))
        isolate_sequence = isolate_records[0].seq
        isolate_gbk_fp = input_dir + '/' + isolate_gbk
        #isolate_gbk_fp = input_dir + isolate_id + '/annomerge/' + isolate_gbk 
        #isolate_gbk_fp = GROUPHOME + '/data/pangenome/gff-copies/all/temp-gbks/1-0006.gbk'
        test_isolate_recs = list(SeqIO.parse(isolate_gbk_fp, 'genbank'))
        output_isolate_recs = [r for r in test_isolate_recs]
        isolate_features = test_isolate_recs[0].features[:]
        output_isolate_recs[0].features = []
        prev_feature_list = []
        num_overlaps = 0
        positions_to_be_resolved = []
        resolve_pairs = []
        for feature in isolate_features:
            if feature.type != 'CDS':
                continue
            if len(prev_feature_list) == 0:
                prev_feature_list = [int(feature.location.start), int(feature.location.end), int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                prev_feature = feature
                continue
            if int(feature.location.start) <= prev_feature_list[1] and int(feature.location.start) >= prev_feature_list[0] and \
                (int(feature.location.start) == prev_feature_list[0] or int(feature.location.end) == prev_feature_list[1]):
                #print(prev_feature_list)
                #print(int(feature.location.start))
                #print('Prev Feature: ' + prev_feature_list[3])
                #print('Current Feature: ' + str(feature.qualifiers['gene'][0]))
                num_overlaps += 1
                positions_to_be_resolved.append((prev_feature_list[0], prev_feature_list[1], prev_feature_list[2]))
                prev_feature_list = [int(feature.location.start), int(feature.location.end), int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                positions_to_be_resolved.append((int(feature.location.start), int(feature.location.end), int(feature.location.strand)))
                resolve_pairs.append([prev_feature, feature])
                prev_feature = feature
            elif int(feature.location.end) <= prev_feature_list[1] and int(feature.location.end) >= prev_feature_list[0] and \
                (int(feature.location.start) == prev_feature_list[0] or int(feature.location.end) == prev_feature_list[1]):
                #print(prev_feature_list)
                #print(int(feature.location.end))
                #print('Prev Feature: ' + prev_feature_list[3])
                #print('Current Feature: ' + str(feature.qualifiers['gene'][0]))
                num_overlaps += 1
                positions_to_be_resolved.append((prev_feature_list[0], prev_feature_list[1], prev_feature_list[2]))
                prev_feature_list = [int(feature.location.start), int(feature.location.end), int(feature.location.strand), str(feature.qualifiers['gene'][0])]
                positions_to_be_resolved.append((int(feature.location.start), int(feature.location.end), int(feature.location.strand)))
                resolve_pairs.append([prev_feature, feature])
                prev_feature = feature
            else:
                prev_feature_list = [int(feature.location.start), int(feature.location.end), int(feature.location.strand), str(feature.qualifiers['gene'][0])]
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
            #print('RATT')
            #print(ratt_annotation)
            #print('Prokka')
            #print(prokka_annotation)
            take_ratt = pick_best_hit(isolate_id, ratt_annotation, prokka_annotation)
            #if = ratt_annotation
            #feature_start = int(f.location.start)
            #feature_stop = int(f.location.end)
            #feature_strand = int(f.location.strand)
            #locus = f.qualifiers['locus_tag'][0]
            #if feature_strand == 1:
            #    if feature_start == 0:
            #        prom_end = len(source_seq)
            #        prom_start = len(source_seq) - 40
            #        prom_seq = str(isolate_sequence[prom_start:prom_end])
            #    elif feature_start < 40:
            #        prom_end = feature_start
            #        prev_prom_len = 40 - prom_end
            #        prom_start = len(source_seq) - prev_prom_len
            #        prom_seq = str(isolate_sequence[prom_start:len(source_seq)]) + str(source_seq[0:prom_end])
            #    else:
            #        prom_end = feature_start
            #        prom_start = feature_start - 40
            #        prom_seq = str(isolate_sequence[prom_start:prom_end])
            #else:
            #    prom_start = feature_stop
            #    prom_end = feature_stop + 40
            #    prom_seq = str(isolate_sequence[prom_start:prom_end])
            #try:
            ###    blast_to_rv_prom = NcbiblastnCommandline(subject=h37rv_prom_fp_dict[locus], outfmt='"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps"')
            #    stdout, stderr = blast_to_rv_prom(stdin=str(prom_seq))
            #    prom_blast_elements = stdout.split('\n')
            #    for line in prom_blast_elements:
            #        if line.startswith('#') or len(line) <= 1:
            #            continue
            #        blast_results = line.strip().split('\t')
            #        if float(blast_results[2]) == 100.0 and int(blast_results[3]) == 40:
            #            prom_mutation = False
            #            break
            #        else:
            #            prom_mutation = True
            #            break
            #except KeyError:
            #    print("Can't find gene: " + locus)
            #    print("WARNING: Taking RATT annotation without checking promoter region for the above gene.")
            #    prom_mutation = False
            if take_ratt:
                output_isolate_recs[0].features.append(ratt_annotation)
            else:
                output_isolate_recs[0].features.append(prokka_annotation)
            #print('Prom Mutation')
            #print(take_ratt)
            #print(blast_results)
            #print(h37rv_prom_fp_dict[locus])
            #print(prom_seq)
            #print('RATT')
            #print(ratt_annotation)
            #print('Prokka')
            #print(prokka_annotation)
        #print(isolate_id)
        #print('Num overlaps = ' + str(num_overlaps))
        #print(positions_to_be_resolved)
        #print(len(resolve_pairs))
        #print(len(isolate_features))
        #print(len(output_isolate_recs[0].features))
        ordered_feats = get_ordered_features(output_isolate_recs[0].features)
        output_isolate_recs[0].features = ordered_feats[:]
        #print(len(output_isolate_recs[0].features))
        SeqIO.write(output_isolate_recs, output_fp, 'genbank')
    for temp_file in h37rv_prom_fp_dict.values():
        os.unlink(temp_file)


if __name__ == "__main__":
    main()
