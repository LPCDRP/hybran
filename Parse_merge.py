#!/usr/bin/env python2.7

import sys
import csv
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import collections
from numpy import median
import os
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
unique_east_asian_merged_genes = []
unique_east_african_indian_merged_genes = []
unique_euro_american_merged_genes = []
unique_indo_oceanic_merged_genes = []
unique_lineage7_merged_genes = []

#Creates a list of filepaths for the gbf files to parse
def get_file_path(DIR):
    filelist = []
    for iso_dir in DIR:
        if iso_dir.startswith('1-') or iso_dir.startswith('2-') or iso_dir.startswith('3-') or iso_dir.startswith('4-') or iso_dir.startswith('SEA'):
            gbf = '/grp/valafar/data/depot/annotation/' + iso_dir + '/' + 'annomerge/merged_genes.gbf'
            filelist.append(gbf)
    return filelist
annomerge_filepaths = get_file_path(root_directory)


#creates a dictionary of each isolate as the key and each merged gene it contains as a value
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
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') or line.startswith('SEA'):
            values = line.split("\t")
            east_asian_list.append(values[0])
    if 'East-African' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') or line.startswith('SEA'):
            values = line.split("\t")
            east_african_indian_list.append(values[0])
    if 'Euro' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') or line.startswith('SEA'):
            values = line.split("\t")
            euro_american_list.append(values[0])
    if 'Indo' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') or line.startswith('SEA'):
            values = line.split("\t")
            indo_oceanic_list.append(values[0])
    if 'Lineage' in line:
        if line.startswith('1-') or line.startswith('2-') or line.startswith('3-') or line.startswith('4-') or line.startswith('SEA'):
            values = line.split("\t")
            lineage7_list.append(values[0])

#Creates a list for each lineage that contains the merged_genes specific to that lineage
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
for gene in flattened_east_asian_gene_list:
    if gene not in flattened_east_african_indian_gene_list and gene not in flattened_euro_american_gene_list and gene not in flattened_indo_oceanic_gene_list and gene not in flattened_lineage7_gene_list:
        if gene not in unique_east_asian_merged_genes:
            unique_east_asian_merged_genes.append(gene)
for gene in flattened_east_african_indian_gene_list:
    if gene not in flattened_indo_oceanic_gene_list and gene not in flattened_lineage7_gene_list and gene not in flattened_euro_american_gene_list and gene not in flattened_east_asian_gene_list:
        if gene not in unique_east_african_indian_merged_genes:
            unique_east_african_indian_merged_genes.append(gene)
for gene in flattened_euro_american_gene_list:
    if gene not in flattened_indo_oceanic_gene_list and gene not in flattened_lineage7_gene_list and gene not in flattened_east_african_indian_gene_list and gene not in flattened_east_asian_gene_list:
        if gene not in unique_euro_american_merged_genes:
            unique_euro_american_merged_genes.append(gene)
for gene in flattened_indo_oceanic_gene_list:
    if gene not in flattened_euro_american_gene_list and gene not in flattened_lineage7_gene_list and gene not in flattened_east_african_indian_gene_list and gene not in flattened_east_asian_gene_list:
        if gene not in unique_indo_oceanic_merged_genes:
            unique_indo_oceanic_merged_genes.append(gene)
for gene in flattened_lineage7_gene_list:
    if gene not in flattened_indo_oceanic_gene_list and gene not in flattened_euro_american_gene_list and gene not in flattened_east_african_indian_gene_list and gene not in flattened_east_asian_gene_list:
        if gene not in unique_lineage7_merged_genes:
            unique_lineage7_merged_genes.append(gene)

#Creates list of merged genes present across all lineages
for gene in all_merged_gene_list:
    if gene in flattened_euro_american_gene_list and gene in flattened_indo_oceanic_gene_list and gene in flattened_east_african_indian_gene_list and gene in flattened_east_asian_gene_list:
        if gene not in in_all_lineages_list:
            in_all_lineages_list.append(gene)

# Returns list of lists [merged_gene and proportion of isolates with that merged_gene within specified lineage.
# Arguments: list of all merged_genes in specified lineage, list containing all isolates within specified lineage, list
# of unique merged_genes in specified lineage
def get_lineage_proportion(gene_list,lineage_list,unique_gene_list):
    gene_proportion_list = []
    for gene in gene_list:
        gene_count = gene_list.count(gene)
        gene_proportion = (str(gene_count) + '/' + str(len(lineage_list)))
        if gene in unique_gene_list:
            if gene not in gene_proportion_list:
                gene_proportion_list.append([gene, gene_proportion])
    return gene_proportion_list

euro_proportion = get_lineage_proportion(flattened_euro_american_gene_list,euro_american_list,unique_euro_american_merged_genes)
east_asian_proportion = get_lineage_proportion(flattened_east_asian_gene_list,east_asian_list,unique_east_asian_merged_genes)
east_african_proportion = get_lineage_proportion(flattened_east_african_indian_gene_list,east_african_indian_list,unique_east_african_indian_merged_genes)
indo_oceanic_proportion = get_lineage_proportion(flattened_indo_oceanic_gene_list,indo_oceanic_list,unique_indo_oceanic_merged_genes)
lineage7_proportion = get_lineage_proportion(flattened_lineage7_gene_list,lineage7_list,unique_lineage7_merged_genes)


# Writes files for each lineage containing merged_genes and the proportion of isolates within that lineage containing that merged_gene
"""with open('unique_east_asian.csv','w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([proportion]) for proportion in east_asian_proportion]
with open('unique_east_african_indian.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([proportion]) for proportion in east_african_proportion]
with open('unique_euro_american.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([proportion]) for proportion in euro_proportion]
with open('unique_indo_oceanic.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([proportion]) for proportion in indo_oceanic_proportion]
with open('unique_lineage7.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([proportion]) for proportion in lineage7_proportion]
with open('in_all_lineages.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow([gene]) for gene in in_all_lineages_list]"""