#!/bin/make -f

roary ?= roary
annotation_dir ?= annotation/
NPROC ?= 1

all: test-weights/clustered_proteins

gffs/gff-done:
	mkdir -p gffs && \
	for gff in $$(ls $(annotation_dir)/*/annomerge/*.gff);
	do
		isolate=$$(basename $$gff | cut -f1 -d.);
		awk '/sequence-region L_contig000001/{flag=1;next}/sequence-region L_contig000002/{flag=0}flag' $$gff | \
		grep -v "##gff-version 3" | \
		sed '1i ##gff-version 3' | \
		sed "s/L_contig/$$isolate-L_contig/g" | \
		awk -F'\t' '{OFS="\t"}{if($$3 !~/RNA/){print $$0}}' > $(dir $@)/$$isolate.gff;
	done && touch $@

test-weights/clustered_proteins: gffs/gff-done
	$(roary) -v -z -s -e -r -f $(dir $@) -p $(NPROC) gffs/*.gff

.ONESHELL:
