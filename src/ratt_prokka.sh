#!/bin/sh

references=$1
fasta=$2
isolate=$3
ref_cds=$4
nproc=$5

mkdir -p prokka prokka-noreference ratt annomerge

if [ ! -f ratt/ratt-home ]
then
	printf "cd ratt/\n
	$RATT_HOME/start.ratt.sh $references $fasta $isolate Strain && touch ratt-done\n
	cd ..\n"
	cd ratt/
	$RATT_HOME/start.ratt.sh $references $fasta $isolate Strain && touch ratt-done
	cd ..
	wait
fi

printf "
prokka    --genus Mycobacterium \n
          --kingdom bacteria \n
          --rfam \n
	  --proteins $ref_cds \n
          --rnammer \n
          --gram pos \n
          --usegenus \n
          --cpus $nproc \n
          --outdir prokka \n
          --prefix $isolate \n
          --force \n
          --centre C \n
          --locustag L \n
	  --quiet \n
      $fasta\n
"
prokka    --genus Mycobacterium \
          --kingdom bacteria \
          --rfam \
	  --proteins $ref_cds \
          --rnammer \
          --gram pos \
          --usegenus \
          --cpus $nproc \
          --outdir prokka \
          --prefix $isolate \
          --force \
          --centre C \
          --locustag L \
	  --quiet \
      $fasta
wait

printf "
prokka    --genus Mycobacterium \n
          --kingdom bacteria \n
          --rfam \n
          --rnammer \n
          --gram pos \n
          --usegenus \n
          --cpus $nproc \n
          --outdir prokka \n
          --prefix $isolate \n
          --force \n
          --centre C \n
          --locustag L \n
	  --quiet \n
      $fasta\n
"

prokka    --genus Mycobacterium \
          --kingdom bacteria \
          --rfam \
          --rnammer \
          --gram pos \
          --usegenus \
          --cpus $nproc \
          --outdir prokka \
          --prefix $isolate \
          --force \
          --centre C \
          --locustag L \
	  --quiet \
      $fasta
wait

sed -i 's/ ; ; ; ; ;/; SV 1; circular; genomic DNA; HTG; PRO;/g; s/order(join/join/g; s/complement(order/complement/g' ratt/*.final.embl
annomerge -i $isolate -o annomerge/$isolate.gbk -l annomerge/$isolate.log > annomerge/annomerge.log && seqret $isolate.gbk $isolate.gff -feature -osf gff
wait
cd ..



