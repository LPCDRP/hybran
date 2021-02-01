#!/bin/sh

nproc=$1
emapper_loc=$2
hybran_tmp_dir=$3

emapper.py --override --tax_scope bactNOG --data_dir $emapper_loc -i eggnog_seqs.fasta --output $hybran_tmp_dir/eggnog-mapper-annotations/mtb_diamond -m diamond --go_evidence experimental --seed_ortholog_evalue 0.000001 --report_orthologs --cpu $nproc --temp_dir $hybran_tmp_dir

emapper.py --override --tax_scope bactNOG --data_dir $emapper_loc -i eggnog_seqs.fasta --output $hybran_tmp_dir/eggnog-mapper-annotations/mtb_hmm -m hmmer -d bactNOG --hmm_evalue 0.000001 --go_evidence experimental  --seed_ortholog_evalue 0.000001 --report_orthologs --cpu $nproc --temp_dir $hybran_tmp_dir

