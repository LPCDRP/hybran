#!/usr/bin/env python2.7
# Author: Deepika Gunasekaran
# Title: Update annotation for novel genes from H37Rv
# Description: This program takes as input, the genbank file from annomerge outputs, compiles a list of CDS candidates
# that do not have an annotated gene name and performs a blastp with H37Rv protein fasta file. If the amino acid
# identity is greater than a certain threshold (default: 85%), this script updates the annotation of the genbank file
# with the corresponding H37Rv hit.


