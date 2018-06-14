#!/usr/bin/env python2.7

import fileinput
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


root_directory = os.listdir('/grp/valafar/data/annotation')

#Returns a list of .gbk files
def get_file_path(DIR):
    filelist = []
    for iso_dir in DIR:
        if iso_dir.endswith('.gbk') and not iso_dir.startswith('A') and not iso_dir.startswith('H37'):
            gbf = '/grp/valafar/data/annotation/' + iso_dir
            filelist.append(gbf)
    return filelist

#Returns a list of lists [isolate_id, locus_tag, location start/end]
domain_list = []
gbf_files_needed = get_file_path(root_directory)
for gbf in gbf_files_needed:
    gbf_string = str(gbf).split('/')
    isolate_id = gbf_string[5]
    isolate_id1 = isolate_id[:-4]
    for record in SeqIO.parse(open(gbf, 'r'), "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "inference" in feature.qualifiers:
                inference_value = feature.qualifiers['inference']
                if "Uniprot" not in str(inference_value) and "HAMAP" in str(inference_value):
                    locus_tag = (str(feature.qualifiers['locus_tag'])).replace("[",'').replace("]",'')
                    if 'L' not in locus_tag and 'Rv' not in locus_tag:
                        domain_list.append([isolate_id1, locus_tag, feature.location.start.position, feature.location.end.position])


#Writes list to file, each sublist as row
"""with open("domain.csv",'wb') as domain_file:
    for domain in domain_list:
        f = csv.writer(domain_file)
        f.writerow([domain])"""

