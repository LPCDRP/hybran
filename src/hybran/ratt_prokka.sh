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

mkdir -p prokka ratt annomerge

if [ ! -f ratt/ratt-done ]
then
	cd ratt/;
	ratt -p "$isolate" -t "$ratt_transfer_type" "$references" "$fasta"
	if [ $? -ne 0 ]
	then
	    exit 1
	fi
	touch ratt-done;
	cd ..
fi

wait

if [ ! -f prokka/"$isolate".gbk ]
then
    prokka \
	$prokka_flags \
	--rfam \
	--rnammer \
	--cpus "$nproc" \
	--coverage "$qcov" \
	--gcode "$gcode" \
	--outdir prokka \
	--prefix "$isolate" \
	--force \
	--centre C \
	--locustag L \
	--quiet \
	"$fasta";
    if [ $? -ne 0 ]
    then
	exit 2
    fi
fi

wait
