#!/bin/sh

nproc=$1
emapper_loc=$2
hybran_tmp_dir=$3
ref_taxon="$4"

emapper.py \
    --override \
    --target_taxa "$ref_taxon"  \
    --data_dir $emapper_loc \
    -i $hybran_tmp_dir/eggnog_seqs.fasta \
    --output $hybran_tmp_dir/eggnog-mapper-annotations/generics_diamond \
    -m diamond \
    --go_evidence experimental \
    --seed_ortholog_evalue 0.000001 \
    --report_orthologs \
    --temp_dir "$hybran_tmp_dir" \
    --cpu $nproc

emapper.py \
    --override \
    --target_taxa "$ref_taxon" \
    --data_dir $emapper_loc \
    -i $hybran_tmp_dir/eggnog_seqs.fasta \
    --output $hybran_tmp_dir/eggnog-mapper-annotations/generics_hmm \
    -m hmmer \
    -d 2 \
    --evalue 0.000001 \
    --go_evidence experimental \
    --seed_ortholog_evalue 0.000001 \
    --report_orthologs \
    --temp_dir "$hybran_tmp_dir" \
    --cpu $nproc
