import os
import logging
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def grep_seqs(gff):
    gff_lines = []
    cmd = ['grep', 'translation=', gff]
    translations = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in translations.stdout:
        try:
            gff_lines.append(line)
        except TypeError:
            continue
    return gff_lines


def create_fasta(directory):
    """Writes cds_seqs.fasta, which contains the CDSs and protein sequences
    from the reference and input genome(s) GFF files."""
    logger = logging.getLogger('CreateFASTA')
    seqs = []
    gff_gene_dict = {}
    for gff in directory:
        raw_out = grep_seqs(gff)
        gff_name = gff.split('/')[-1].split('.')[0]
        for line in raw_out:
            gff_id = gff_name + '-' + [j.split('=')[1] for i in line.split('\t')
                                       if i.startswith('ID=') for j in i.split(';')][0]
            gene = None
            for i in line.split(';'):
                if i.startswith('gene='):
                    gene = [i.split('=')[1].rstrip('\n')][0]
                if i.startswith('locus_tag='):
                    locus_tag = [i.split('=')[1].rstrip('\n')][0]
            if not gene:
                gene = locus_tag
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
