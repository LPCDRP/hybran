#!/bin/make -f

IN_GENOME ?= isolate.fasta
H37Rv_FASTA ?= H37Rv-CDS.fasta
NPROC ?= 1

export RATT_CONFIG := $(RATT_HOME)/RATT.config
H37Rv_EMBL_PATH := $(RATT_HOME)/embl/
ratt = $(RATT_HOME)/start.ratt.sh

all: annomerge/%.gbk

%/annomerge/%.gbk: %/ratt/ratt-done %/prokka/%.gbf %/prokka-noreference/%.gbf
	mkdir -p $*/annomerge
	cd $*/annomerge
	sed -i 's/ ; ; ; ; ;/; SV 1; circular; genomic DNA; HTG; PRO;/g; s/order(join/order/g; s/complement(order/complement/g' ../ratt/*.final.embl
	annomerge -i $* -o $*.gbk -l $*.log > annomerge.log \
	 && seqret $*.gbk $*.gff -feature -osf gff

%/ratt/ratt-done: %.fasta %/ $(H37Rv_EMBL_PATH)
	mkdir -p ratt
	cd ratt
	$(ratt) $(word 2, $^) $< $* Strain && touch ratt-done

%/prokka/%.gbf: %.fasta %/
	prokka    --genus Mycobacterium \
	          --kingdom bacteria \
	          --rfam \
		  --proteins $(H37Rv_FASTA) \
	          --rnammer \
	          --gram pos \
	          --usegenus \
	          --cpus $(NPROC) \
	          --outdir $(dir $@) \
	          --prefix $* \
	          --force \
	          --centre C \
	          --locustag $* \
	      $<
	touch $@

%/prokka-noreference/%.gbf: %.fasta %/
	prokka    --genus Mycobacterium \
	          --kingdom bacteria \
	          --rfam \
	          --rnammer \
	          --gram pos \
	          --usegenus \
	          --cpus $(NPROC) \
	          --outdir $(dir $@) \
	          --prefix $* \
	          --force \
	          --centre C \
	          --locustag $* \
	      $<
	touch $@

%/:
	mkdir -p $@

.DELETE_ON_ERROR:
.SECONDARY:
.ONESHELL:
