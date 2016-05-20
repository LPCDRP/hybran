#!/usr/bin/make -f

SAMPLE ?= A6

basepath ?= $(GROUPHOME)/data/depot/annotation
VPATH += $(addprefix $(basepath)/, $(SAMPLE) $(SAMPLE)/prokka $(SAMPLE)/Reference)

# At least the prokkAnn won't be consistently named across samples...
prokkAnn ?= L_12152015.gff
rattAnn ?= $(SAMPLE).fasta.1_quiver_quiver_quiver.final.embl
rattNonsyn ?= $(SAMPLE).fasta.H37Rv_NC000962_3.Mutations.gff

.PHONY: all
all: $(SAMPLE).merged.embl

%.merged.embl: %.intersect.bed $(rattAnn).in
	python annomerge.py -b $(word 1,$^) -e $(word 2,$^) -o $@

%.intersect.bed: $(prokkAnn).in $(rattNonsyn).in
	bedtools intersect -a $(word 1,$^) -b $(word 2,$^) -f 1 > $@

# Sam's code to convert RATT's output into valid EMBL
%.embl.in: %.embl
	sed  's/ ; ; ; ; ;/; SV 1; circular; genomic DNA; HTG; PRO;/g' $< \
	| sed 's/order(join/order/g' \
	| sed 's/complement(order/complement/g' \
	> $@

%.Mutations.gff.in: %.Mutations.gff
	python formatRattForBedtools.py $< $@

%.gff.in: %.gff
	python formatProkkaForBedtools.py $< $@
