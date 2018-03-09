#Makefile to organize creation of reference genome embl file for annomerge
#Uses tbl file and csv files with TSS data
#call like this:
#make -i -f tbl2tssembl.mk
.ONE_SHELL:
#Variables
finalembl ?= ittas-man-h37-tss-merge.embl

#reference feature table
tblsource ?= $(GROUPHOME)/resources/H37Rv-hypotheticome.tbl
#tblsource ?= /home/dconkleg/H37Rv-mannotation-computation.tbl
fastasource := $(GROUPHOME)/resources/H37Rv-NC_000962.3.fasta
templatesource := $(GROUPHOME)/resources/tbl2asn-template.sbt

#path to TSS data
tssdir := $(GROUPHOME)/resources/tss-csv-data/

tempdir = temp-tbl2asn-in/
temptable = sequence.tbl           
tempfasta = sequence.fsa
tempasn = sequence.sqn
tempembl = temp-ittas-man-h37.embl

fastaheader = >NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome
tblheader = >Feature NC_000962.3 Table1


#Recipies
all : $(finalembl)

#adds TSS data to embl file of reference genome
$(finalembl) : $(tempdir)$(tempembl)
	#python tssmerger.py -r '$<' -o '$(finalembl)' -c '$(tssdir)'  #test version
	tssmerger -r '$<' -o '$(finalembl)' -c '$(tssdir)'
	rm -rf $(tempdir)

#create embl from asn 
$(tempdir)$(tempembl) : $(tempdir)$(tempasn)
	asn2gb -i $< -o $@ -f e

#create asn from tbl
$(tempdir)$(tempasn) : $(templatesource) $(tempdir)$(tempfasta) $(tempdir)$(temptable)
	tbl2asn -k m -t $< -p $(tempdir) -r $(tempdir)

#create tbl file from source
$(tempdir)$(temptable) : $(tblsource) $(tempdir)
	#cp $< $(tempdir)tablemid1
	sed 's/EC number/EC_number/' $<	> $(tempdir)tablemid1;\
	sed -i '1d' $(tempdir)tablemid1;\
	sed -i '1s/^/>Feature NC_000962.3 Table1\n/' $(tempdir)tablemid1; \
	mv $(tempdir)tablemid1 $@

#create fasta file from source
$(tempdir)$(tempfasta) : $(fastasource) $(tempdir)
	cp $< $(tempdir)sedin
	sed '1d' $(tempdir)sedin > $(tempdir)sedout
	sed -i '1s/^/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome\n/' $(tempdir)sedout
	rm $(tempdir)sedin
	mv $(tempdir)sedout $@
	
#create temp directory
$(tempdir): 
	mkdir -p $@
