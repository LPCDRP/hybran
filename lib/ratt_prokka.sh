#!/bin/sh

references=$1
fasta=$2
isolate=$3
ref_cds=$4
nproc=$5

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	cd ratt/ && \ 
	$RATT_HOME/start.ratt.sh $references $fasta $isolate Strain && \ 
	touch ratt-done && \ 
	cd ..
fi

wait

if [ ! -f prokka/$isolate.gbk ]
then
	printf "prokka --genus Mycobacterium --kingdom bacteria --rfam --proteins $ref_cds --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka --prefix $isolate --force --centre C --locustag L --quiet $fasta\n"
	prokka --genus Mycobacterium --kingdom bacteria --rfam --proteins $ref_cds --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka --prefix $isolate --force --centre C --locustag L --quiet $fasta
else
	printf "prokka/$isolate.gbk already exists\n"
fi

wait

if [ ! -f prokka-noreference/$isolate.gbk ]
then
	printf "prokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka-noreference --prefix $isolate --force --centre C --locustag L --quiet $fasta\n"

prokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus $nproc --outdir prokka-noreference --prefix $isolate --force --centre C --locustag L --quiet $fasta
else
	printf "prokka-noreference/$isolate.gbk already exists\n"
fi

wait

if [ ! -f annomerge/$isolate.gbk ]
then
	sed -i 's/ ; ; ; ; ;/; SV 1; circular; genomic DNA; HTG; PRO;/g; s/order(join/join/g; s/complement(order/complement/g' ratt/*.final.embl
	annomerge -i $isolate -o annomerge/$isolate.gbk -l annomerge/$isolate.log > annomerge/annomerge.log && seqret annomerge/$isolate.gbk annomerge/$isolate.gff -feature -osf gff
else
	printf "annomerge/$isolate.gbk already exists\n"
fi

wait
cd ..

