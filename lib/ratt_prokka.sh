#!/bin/sh

references=$1
fasta=$2
isolate=$3
ref_cds=$4
nproc=$5

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	printf "\t\t$RATT_HOME/start.ratt.sh $references $fasta $isolate Strain\n"
	cd ratt/;
	$RATT_HOME/start.ratt.sh $references $fasta $isolate Strain; 
	touch ratt-done;
	cd ..
else
	printf "Annotations already transferred to $isolate\n"
fi

wait

if [ ! -f prokka/$isolate.gbk ]
then
	printf "\t\tprokka --genus Mycobacterium --kingdom bacteria --rfam --proteins $ref_cds --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka --prefix $isolate --force --centre C --locustag L --quiet $fasta\n"
	prokka --genus Mycobacterium --kingdom bacteria --rfam --proteins $ref_cds --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka --prefix $isolate --force --centre C --locustag L --quiet $fasta
else
	printf "\t\tprokka/$isolate.gbk already exists\n"
fi

wait

if [ ! -f prokka-noreference/$isolate.gbk ]
then
	printf "\t\tprokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka-noreference --prefix $isolate --force --centre C --locustag L --quiet $fasta\n"

	prokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka-noreference --prefix $isolate --force --centre C --locustag L --quiet $fasta
else
	printf "\t\tprokka-noreference/$isolate.gbk already exists\n"
fi

wait

if [ ! -f $isolate.gbk ]
then
	printf "\t\tannomerge -i $isolate -o $isolate.gbk\n"
	printf "lib/annomerge.py -i $isolate -o $isolate.gbk -l annomerge/$isolate.log -fp ./ -prot $ref_cds > annomerge/annomerge.log\n"
	sed -i 's/ ; ; ; ; ;/; SV 1; circular; genomic DNA; HTG; PRO;/g; s/order(join/join/g; s/complement(order/complement/g' ratt/*.final.embl
	python lib/annomerge.py -i $isolate -o $isolate.gbk -l $isolate/annomerge/$isolate.log -fp ./ -prot $ref_cds > $isolate/annomerge/annomerge.log
else
	printf "\t\t$isolate.gbk already exists\n"
fi

if [ ! -f $isolate.gff ]
then
	printf "\t\tseqret $isolate.gbk $isolate.gff -feature -osf gff\n"
	seqret $isolate.gbk $isolate.gff -feature -osf gff 2> /dev/null
else
	printf "\t\t$isolate.gff already exists\n"
fi

wait
cd ..

