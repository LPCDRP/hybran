#!/bin/sh

isolate_fasta=$1

mkdir -p prodigal

cd prodigal/ && \
prodigal -i $isolate_fasta -o reference_prodigal -a reference_prodigal_proteome.faa && \
cd ..


