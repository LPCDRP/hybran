#!/bin/sh

references=$1
fasta=$2
isolate=$3
ref_cds=$4
nproc=$5
qcov=$6
gcode=$7
ratt_transfer_type="$8"
prokka_flags="$9"

set -x

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	cd ratt/;
	ratt -p "$isolate" -t "$ratt_transfer_type" "$references" "$fasta"
	touch ratt-done;
	cd ..
fi

wait

if [ ! -f prokka/"$isolate".gbk ]
then
    prokka \
	$prokka_flags \
	--rfam \
	--proteins "$ref_cds" \
	--rnammer  \
	`# rawproduct: prevents removal of gene name assignments when reference product field isn't specific` \
	--rawproduct \
	--coverage 0 \
	--gcode "$gcode" \
	--cpus "$nproc" \
	--outdir prokka \
	--prefix "$isolate" \
	--force \
	--centre C \
	--locustag L \
	--quiet \
	"$fasta";
fi

wait

if [ ! -f prokka-noreference/"$isolate".gbk ]
then
    prokka \
	$prokka_flags \
	--rfam \
	--rnammer \
	--cpus "$nproc" \
	--coverage "$qcov" \
	--gcode "$gcode" \
	--outdir prokka-noreference \
	--prefix "$isolate" \
	--force \
	--centre C \
	--locustag L \
	--quiet \
	"$fasta";
fi

wait
