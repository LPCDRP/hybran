#!/bin/sh

references=$1
fasta=$2
isolate=$3
ref_cds=$4
nproc=$5

set -x

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	cd ratt/;
	ratt -p "$isolate" -t Strain "$references" "$fasta"
	touch ratt-done;
	cd ..
fi

wait

if [ ! -f prokka/"$isolate".gbk ]
then
	prokka --genus Mycobacterium --kingdom bacteria --rfam --proteins "$ref_cds" --rnammer --gram pos --usegenus \
	--cpus "$nproc" --outdir prokka --prefix "$isolate" --force --centre C --locustag L --quiet "$fasta";
fi

wait

if [ ! -f prokka-noreference/"$isolate".gbk ]
then
	prokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus "$nproc" \
	--outdir prokka-noreference --prefix "$isolate" --force --centre C --locustag L --quiet "$fasta";
fi

wait
