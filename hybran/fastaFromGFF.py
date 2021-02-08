import os
import logging
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def grep_seqs(gff):
    cmd = ['grep', 'translation=', gff]
    translations = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return [line for line in translations.stdout]


def create_fasta(directory):
    logger = logging.getLogger('CreateFASTA')
    seqs = []
    gff_gene_dict = {}
    for gff in directory:
        raw_out = grep_seqs(gff)
        gff_name = gff.split('/')[-1].split('.')[0]
        for line in raw_out:
            gff_id = gff_name + '-' + [j.split('=')[1] for i in line.split('\t')
                                       if i.startswith('ID=') for j in i.split(';')][0]
            gene = [i.split('=')[1].rstrip('\n') for i in line.split(';') if i.startswith('gene=')][0]
            gff_gene_dict[gff_id] = gene
            translation = [i.split('=')[1] for i in line.split(';') if i.startswith('translation=')][0]
            record = SeqRecord(Seq(translation.rstrip('\n')),
                               id=gff_id,
                               description='')
            seqs.append(record)
    logger.info('Writing FASTA to cds_seqs.fasta')
    with open('cds_seqs.fasta', 'w') as out:
        for s in seqs:
            SeqIO.write(s, out, 'fasta')
    return gff_gene_dict
