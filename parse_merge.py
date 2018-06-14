#!/usr/bin/env python2.7

import sys
import csv
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import collections
from numpy import median
import os, os.path
from Bio import GenBank
import fileinput

all_merged_gene_list = []
in_all_lineages_list = []
isolate_id_list = []
gene_dict = {}
root_directory = os.listdir('/grp/valafar/data/depot/annotation')
east_african_indian_list = []
east_asian_list = []
euro_american_list = []
indo_oceanic_list = []
lineage7_list = []
east_african_indian_gene_dict = {}
east_asian_gene_dict = {}
euro_american_gene_dict = {}
indo_oceanic_gene_dict = {}
lineage7_gene_dict = {}
east_african_indian_gene_list = []
east_asian_gene_list = []
euro_american_gene_list = []
indo_oceanic_gene_list = []
lineage7_gene_list = []
flattened_east_african_indian_gene_list = []
flattened_east_asian_gene_list = []
flattened_euro_american_gene_list = []
flattened_indo_oceanic_gene_list = []
flattened_lineage7_gene_list = []




#Creates a list of filepaths for the gbf files to parse
def get_file_path(DIR):
    filelist = []
    for iso_dir in DIR:
        if iso_dir.startswith('1-') or iso_dir.startswith('2-') or iso_dir.startswith('3-') or iso_dir.startswith('4-') \
                or iso_dir.startswith('SEA'):
            if os.path.exists('/grp/valafar/data/depot/annotation/' + iso_dir + '/' + 'annomerge/merged_genes.gbf'):
                gbf = '/grp/valafar/data/depot/annotation/' + iso_dir + '/' + 'annomerge/merged_genes.gbf'
                filelist.append(gbf)
    return filelist
annomerge_filepaths = get_file_path(root_directory)


#creates a dictionary of each isolate as the key and each merged gene it contains as a value
# Also creates a list of each merged gene
for gbf in annomerge_filepaths:
    parser = GenBank.RecordParser()
    record = parser.parse(open(gbf))
    features = record.features
    gbf_split = gbf.split('/')
    isolate_id = gbf_split[6]
    gene_list = []
    for feature in features:
        gene_comment_list = []
        qualifiers = feature.qualifiers
        for qualifier in qualifiers:
            qvalue = qualifier.value
            if 'merged' in qvalue:
                gene_comment_list.append(qvalue)
                for gene_comment in gene_comment_list:
                    gene_split_comment = gene_comment.split()
                    gene_list.append(gene_split_comment[2])
                    gene_dict[isolate_id] = gene_list
                    all_merged_gene_list.append(gene_split_comment[2])



#Creating lists of every isolate for each lineage
with open('/grp/valafar/metadata/lineage.txt') as lineage_file:
    content = lineage_file.readlines()
for line in content:
    if 'East-Asian' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') \
                or line.startswith('SEA'):
            values = line.split("\t")
            east_asian_list.append(values[0])
    if 'East-African' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') \
                or line.startswith('SEA'):
            values = line.split("\t")
            east_african_indian_list.append(values[0])
    if 'Euro' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') \
                or line.startswith('SEA'):
            values = line.split("\t")
            euro_american_list.append(values[0])
    if 'Indo' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') \
                or line.startswith('SEA'):
            values = line.split("\t")
            indo_oceanic_list.append(values[0])
    if 'Lineage' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') \
                or line.startswith('SEA'):
            values = line.split("\t")
            lineage7_list.append(values[0])


#Creates a dictionary for each lineage that contains the isolates and merged_genes present in that lineage
for isolate, gene in gene_dict.items():
    if isolate in east_asian_list:
        east_asian_gene_dict.update({isolate:gene})
    if isolate in east_african_indian_list:
        east_african_indian_gene_dict.update({isolate:gene})
    if isolate in euro_american_list:
        euro_american_gene_dict.update({isolate:gene})
    if isolate in indo_oceanic_list:
        indo_oceanic_gene_dict.update({isolate:gene})
    if isolate in lineage7_list:
        lineage7_gene_dict.update({isolate:gene})


# Creates a list for each lineage that contains the isolates present in that lineage
# Used for comparison in the get_unique_genes function
for isolate, gene in gene_dict.items():
    if isolate in east_asian_list:
        east_asian_gene_list.append(gene)
    if isolate in east_african_indian_list:
        east_african_indian_gene_list.append(gene)
    if isolate in euro_american_list:
        euro_american_gene_list.append(gene)
    if isolate in indo_oceanic_list:
        indo_oceanic_gene_list.append(gene)
    if isolate in lineage7_list:
        lineage7_gene_list.append(gene)


#Flattens the lists
for sublist in east_asian_gene_list:
    for gene in sublist:
        flattened_east_asian_gene_list.append(gene)
for sublist in east_african_indian_gene_list:
    for gene in sublist:
        flattened_east_african_indian_gene_list.append(gene)
for sublist in euro_american_gene_list:
    for gene in sublist:
        flattened_euro_american_gene_list.append(gene)
for sublist in indo_oceanic_gene_list:
    for gene in sublist:
        flattened_indo_oceanic_gene_list.append(gene)
for sublist in lineage7_gene_list:
    for gene in sublist:
        flattened_lineage7_gene_list.append(gene)


#Creates a list for each lineage that contains merged genes not found in any other lineage
def get_unique_genes(gene_dict,compared_list1,compared_list2,compared_list3,compared_list4):
    unique_dict = {}
    unique_gene_list = []
    for isolate,genes in gene_dict.items():
        for gene in genes:
            if gene not in compared_list1 and gene not in compared_list2 and gene not in compared_list3 and \
                    gene not in compared_list4:
                if gene not in unique_gene_list:
                    unique_dict.setdefault(isolate,[]).append(gene)
                    unique_gene_list.append(gene)

    return unique_dict



unique_east_african_indian_merged_genes = get_unique_genes(east_african_indian_gene_dict,flattened_east_asian_gene_list,
                        flattened_euro_american_gene_list,flattened_indo_oceanic_gene_list,flattened_lineage7_gene_list)
unique_east_asian_merged_genes = get_unique_genes(east_asian_gene_dict,flattened_east_african_indian_gene_list,
                        flattened_euro_american_gene_list,flattened_indo_oceanic_gene_list,flattened_lineage7_gene_list)
unique_euro_american_merged_genes = get_unique_genes(euro_american_gene_dict,flattened_east_african_indian_gene_list,
                        flattened_east_asian_gene_list,flattened_indo_oceanic_gene_list,flattened_lineage7_gene_list)
unique_indo_oceanic_merged_genes = get_unique_genes(indo_oceanic_gene_dict,flattened_east_african_indian_gene_list,
                        flattened_east_asian_gene_list,flattened_euro_american_gene_list,flattened_lineage7_gene_list)
unique_lineage7_merged_genes = get_unique_genes(lineage7_gene_dict,flattened_east_african_indian_gene_list,
                        flattened_east_asian_gene_list,flattened_euro_american_gene_list,flattened_indo_oceanic_gene_list)




def get_counts(unique_gene_dict,all_genes_list):
    count_dict = {}
    unique_gene_list = []
    for isolate,genes in unique_gene_dict.items():
        for gene in genes:
            if gene not in unique_gene_list:
                count_dict.setdefault(isolate,[]).append([gene,all_genes_list.count(gene)])
                unique_gene_list.append(gene)
    return count_dict


ALL_counts = get_counts(gene_dict,all_merged_gene_list)
IO_counts = get_counts(unique_indo_oceanic_merged_genes,all_merged_gene_list)
EUR_counts = get_counts(unique_euro_american_merged_genes,all_merged_gene_list)
EAI_counts = get_counts(unique_east_african_indian_merged_genes,all_merged_gene_list)
EAS_counts = get_counts(unique_east_asian_merged_genes,all_merged_gene_list)
LIN7_counts = get_counts(unique_lineage7_merged_genes,all_merged_gene_list)









# Writes files for each lineage containing merged_genes and the proportion of isolates within that lineage containing that merged_gene
with open('merged_genes.csv','w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Regardless of Lineage'])
    for isolate,genes in ALL_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])
with open('lineage_specific_merged_genes.csv','w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['East African Indian Lineage'])
    for isolate,genes in EAI_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])
    writer.writerow([])
    writer.writerow(['East Asian Lineage'])
    for isolate,genes in EAS_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])
    writer.writerow([])
    writer.writerow(['Euro American Lineage'])
    for isolate,genes in EUR_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])
    writer.writerow([])
    writer.writerow(['Indo Oceanic Lineage'])
    for isolate,genes in IO_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])
    writer.writerow([])
    writer.writerow(['Lineage 7'])
    for isolate,genes in LIN7_counts.items():
        writer.writerow([])
        writer.writerow([isolate])
        for gene in genes:
            writer.writerows([gene])




