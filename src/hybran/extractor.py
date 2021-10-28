import logging
import os
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fastaFromGbk(genbank, out_cds, out_genome):
    """
    Extracts amino acid CDS sequences from a Genbank annotation file

    :param genbank:
    :param out_cds: file name or handle in which to write CDS sequences
    :param out_genome: file name or handle in which to write genome sequence
    :return:
    """
    logger = logging.getLogger('FastaFromGbk')
    seqs = []
    for record in SeqIO.parse(genbank, 'genbank'):
        seq = SeqRecord(record.seq, id=os.path.splitext(os.path.basename(genbank))[0],
                        description='')
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS' and 'pseudogene' not in feature.qualifiers:
                    try:
                        gene = feature.qualifiers['gene'][0]
                    except KeyError:
                        gene = ''
                    try:
                        locus_tag = feature.qualifiers['locus_tag'][0]
                    except KeyError:
                        locus_tag = ''
                    if not locus_tag and not gene:
                        logger.error('No locus tag or gene name for the following Genbank entry. Please check the '
                                     'Genbank file\n' + str(feature))
                        exit(KeyError)
                    seq_record_id = locus_tag + ':' + gene
                    seq_record = SeqRecord(Seq(feature.qualifiers['translation'][0]),
                                           id=seq_record_id,
                                           description='')
                    seqs.append(seq_record)
    SeqIO.write(seqs, out_cds, 'fasta')
    SeqIO.write(seq, out_genome, 'fasta')
    return


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


def fastaFromGffList(gffs, out_cds):
    """
    Extracts amino acid CDS sequences from GFF annotation files

    :param gffs: list of GFF annotation file names
    :param out_cds: file name or handle in which to write CDS sequences
    :return gff_gene_dict: dictionary of CDS sequences by sample-recordID
    """
    logger = logging.getLogger('FastaFromGff')
    seqs = []
    gff_gene_dict = {}
    for gff in gffs:
        raw_out = grep_seqs(gff)
        gff_name = os.path.splitext(os.path.basename(gff))[0]
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

    SeqIO.write(seqs, out_cds, 'fasta')
    return gff_gene_dict
