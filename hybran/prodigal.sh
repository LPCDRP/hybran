#!/bin/sh

isolate_fasta=$1
hybran_tmp_dir=$2

mkdir -p $hybran_tmp_dir/prodigal-test
prodigal -i $isolate_fasta -o $hybran_tmp_dir/prodigal-test/reference_prodigal -a reference_prodigal_proteome.faa
