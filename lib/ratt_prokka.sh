#!/bin/sh

references=$1
prokka_noref_ltag=$2
fasta=$3
isolate=$4
ref_cds=$5
nproc=$6

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	cd ratt/;
	$RATT_HOME/start.ratt.sh "$references" "$fasta" "$isolate" Strain || exit 1
	touch ratt-done;
	cd ..
fi

wait

if [ ! -f prokka/"$isolate".gbk ]
then
  prokka --genus Mycobacterium --kingdom bacteria --rfam --proteins "$ref_cds" --rnammer --gram pos --usegenus --cpus "$nproc" --outdir prokka --prefix "$isolate" --force --centre C --locustag "$isolate" --quiet "$fasta"
fi

wait

if [ ! -f prokka-noreference/"$isolate".gbk ]
then
	prokka --genus Mycobacterium --kingdom bacteria --rfam --rnammer --gram pos --usegenus --cpus "$nproc" --outdir prokka-noreference --prefix "$isolate" --force --centre C --locustag L2 --quiet "$fasta"
fi

wait
