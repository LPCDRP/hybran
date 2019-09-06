#!/bin/sh

isolate_id=$1
isolate_fasta=$2

mkdir -p prodigal

cd prodigal/ && \
prodigal -i $isolate_fasta -o $isolate -a $isolate.faa && \
cd ..


