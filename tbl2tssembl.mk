#Makefile to organize creation of reference genome embl file for annomerge
#Uses tbl file and csv files with TSS data
#call like this:
#make -i -f tbl2tssembl.mk

finalembl ?= ittas-man-h37-tss-merge.embl

#reference feature table
#note: this default tbl filename will need to be updated 
#once Allyssa has finished the latest version
tblsource ?= $(GROUPHOME)/resources/H37Rv-mannotation-computation.tbl      

fastasource := $(GROUPHOME)/resources/H37Rv-NC_000962.3.fasta

#path to TSS data
tssdir := $(GROUPHOME)/resources/tss-csv-data/

templatesource := $(GROUPHOME)/resources/tbl2asn-template.sbt

tempdir = tbl2asn-space/

#these two must have the same prefix
temptable = sequence.tbl           
tempfasta = sequence.fsa

tempasn = sequence.sqn
tempembl = temp-ittas-man-h37.embl

fastaheader = >NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome

all : $(finalembl)

#adds TSS data to embl file of reference genome
$(finalembl) : $(tempdir)$(tempembl)
	tssmerger -r '$<' -o '$(finalembl)' -c '$(tssdir)'
	#python tssmerger.py -r '$<' -o '$(finalembl)' -c '$(tssdir)'  #test version
	rm -rf $(tempdir)

#create embl from asn 
$(tempdir)$(tempembl) : $(tempdir)$(tempasn)
	asn2gb -i $< -o $@ -f e

#create asn from tbl
$(tempdir)$(tempasn) : $(templatesource) $(tempdir)$(tempfasta) $(tempdir)$(temptable)
	tbl2asn -t $< -p $(tempdir) -r $(tempdir)

#create tbl file from source
$(tempdir)$(temptable) : $(tblsource) $(tempdir)
	cp $< $@

#create fasta file from source
$(tempdir)$(tempfasta) : $(fastasource) $(tempdir) 
	cp $< $(tempdir)sedin
	sed '1d' $(tempdir)sedin > $(tempdir)sedout
	sed -i '1s/^/$(fastaheader)\n/' $(tempdir)sedout
	mv $(tempdir)sedout $@
	rm $(tempdir)sedin

#create temp directory
$(tempdir) : 
	mkdir $@
