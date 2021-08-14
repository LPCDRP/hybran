#!/bin/sh

isolate_fasta=$1
outfile=$2

mkdir -p $(dirname "$outfile")
prodigal -i $isolate_fasta -o "$outfile" -a reference_prodigal_proteome.faa
