#!/usr/bin/env python2.7

import os
from Bio import SeqIO

GROUPHOME = os.environ['GROUPHOME']


def check_for_remainder(embl_file):
    locus_tag_list = []
    features = SeqIO.read(embl_file, 'embl').features
    for feature in features:
        if feature.type == 'CDS':
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            if (end - start + 1) % 3 != 0:
                locus_tag_list.append(feature.qualifiers['locus_tag'][0])
    return locus_tag_list
