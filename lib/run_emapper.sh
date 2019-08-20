#!/bin/sh

nproc=$1
emapper_loc=$2

emapper.py 				\ 
	--tax_scope bactNOG		\ 
	--data_dir $emapper_loc 	\ 
	-i eggnog_seqs.fasta 		\ 
	--output mtb_diamond 		\ 
	-m diamond 			\ 
	--go_evidence experimental 	\ 
	--seed_ortholog_evalue 0.000001 \ 
	--report_orthologs 		\ 
	--cpu $nproc

emapper.py 				\ 
	--tax_scope bactNOG		\ 
	--data_dir $emapper_loc 	\ 
	-i eggnog_seqs.fasta 		\ 
	--output mtb_hmm 		\ 
	-m hmmer 			\ 
	-d bactNOG 			\ 
	--hmm_evalue 0.000001 		\ 
	--go_evidence experimental 	\ 
	--seed_ortholog_evalue 0.000001 \ 
	--report_orthologs 		\ 
	--cpu $nproc

