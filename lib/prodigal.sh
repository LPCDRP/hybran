#!/bin/sh

isolate_fasta=$1

mkdir prodigal-test
prodigal -i $isolate_fasta -o prodigal-test/reference_prodigal -a reference_prodigal_proteome.faa
