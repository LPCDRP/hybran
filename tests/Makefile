.SECONDARY:
.DELETE_ON_ERROR:
NPROC ?= 8
RESDIR ?= ${GROUPHOME}/data/depot/hybran/resources
OUTDIR ?= $(shell date +%m_%d)
HYBRANFLAGS ?=

ifeq ($(SUBMIT), 1)
#QSUB ?= qsub -N crazy_ivan -cwd -pe smp $(NPROC) -S /bin/bash -b yes -shell no -V -j yes -sync yes -o $(@D)/run.log
QSUB ?= sbatch -J crazy_ivan -c $(NPROC) --wait -o $(@D)/run.log
else
log = 2>&1 | tee -a $(@D)/run.log
endif

cases = \
 1-0006 \
 1-0006-multiref \
 1-0006-filterratt \
 1-0006-ncbi \
 SRR23086714_FLYE \
 mpox \
 PAK \
 AZ20 \
 AZ20-reannotatedref \


crazy-ivan: $(cases:%=$(OUTDIR)/%/coords.txt) ;
$(cases): %: $(OUTDIR)/%/coords.txt


%/coords.txt: %/final.gbk
	./gbk2coords $< > $@


$(OUTDIR)/1-0006%: genome=$(RESDIR)/1-0006.fasta
$(OUTDIR)/1-0006/final.gbk $(OUTDIR)/1-0006-filterratt/final.gbk $(OUTDIR)/SRR23086714_FLYE/final.gbk: ref=$(RESDIR)/H37Rv-3.3.3.gbk
$(OUTDIR)/1-0006-ncbi/final.gbk: ref=$(RESDIR)/h37rv_ncbi_2023_06_07.gbk
$(OUTDIR)/1-0006-multiref/final.gbk: ref=$(RESDIR)/tb-multiref
$(OUTDIR)/1-0006-filterratt/final.gbk: local_HYBRANFLAGS += --filter-ratt
$(OUTDIR)/SRR23086714_FLYE/final.gbk: genome=$(RESDIR)/SRR23086714_FLYE.fasta
$(OUTDIR)/mpox/final.gbk: genome=$(RESDIR)/UnicyclerMpox.fasta
$(OUTDIR)/mpox/final.gbk: ref=$(RESDIR)/mpox.gbk
$(OUTDIR)/PAK/final.gbk: genome=$(RESDIR)/PAK.fasta
$(OUTDIR)/PAK/final.gbk: ref=/grp/valafar/data/depot/hybran.pub/references/PAO1-refseq-2019.gbk
$(OUTDIR)/AZ20/final.gbk $(OUTDIR)/AZ20-reannotatedref/final.gbk: genome=$(RESDIR)/AZ20.fasta
$(OUTDIR)/AZ20/final.gbk: ref=$(RESDIR)/nissle-hybrid.gbk
$(OUTDIR)/AZ20-reannotatedref/final.gbk: ref=$(RESDIR)/ECOLIN.gbk
%/final.gbk:
	mkdir -p $(@D); $(QSUB) \
	hybran \
	-g $(genome) \
	-r $(ref) \
	--nproc $(NPROC) \
	--verbose \
	$(local_HYBRANFLAGS) $(HYBRANFLAGS) \
	-o $(@D) \
	$(log)
	ln -sf $$(basename $(genome) .fasta).gbk $@
