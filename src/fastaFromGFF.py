import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def grep_seqs(gff):
    cmd = ['grep', 'translation=', gff]
    translations = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return [line for line in translations.stdout]


def create_fasta(directory):
    seqs = []
    gff_gene_dict = {}
    for gff in os.listdir(directory):
        if gff.endswith('.gff'):
            raw_out = grep_seqs(directory + gff)
            for line in raw_out:
                gff_id = [i.split('=')[1] for i in line if i.startswith('ID=')][0]
                gene = [i.split('=')[1] for i in line if i.startswith('gene=')][0]
                gff_gene_dict[gff_id] = gene
                translation = [i.split('=')[1] for i in line if i.startswith('translation=')][0]
                record = SeqRecord(Seq(translation),
                                   id=gff_id,
                                   description='')
                seqs.append(record)
    with open('cds_seqs.fasta', 'w') as out:
        for s in seqs:
            SeqIO.write(s, out, 'fasta')
    return gff_gene_dict
