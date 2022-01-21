import logging
import os
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord


def get_ltag(feature):
    try:
        ltag = feature.qualifiers['locus_tag'][0]
    except KeyError:
        exit("No locus tag for " + str(feature))
    return ltag

def get_gene(feature):
    try:
        gene = feature.qualifiers['gene'][0]
    except KeyError:
        gene = get_ltag(feature)
    return gene

def prokka_faa(feature):
    """
    Create fasta descriptions in Prokka format.
    See https://github.com/tseemann/prokka#fasta-database-format
    :param feature: sequence feature
    :return: string representing fasta record description in prokka format
    """
    # Don't include EC for now.
    # Not sure how to handle multiple EC numbers in Prokka faa format
    # ec = feature.qualifiers['EC_number'][0]
    ec = ''
    try:
        product = feature.qualifiers['product'][0]
    except KeyError:
        product = ''
    return '~~~'.join([ec,get_gene(feature),product])

def fastaFromGbk(genbank, out_cds, out_genome,
                 identify = lambda f: ':'.join([get_ltag(f),
                                                get_gene(f)]),
                 describe = lambda f: '',
                 ):
    """
    Extracts amino acid CDS sequences from a Genbank annotation file

    :param genbank:
    :param out_cds: file name or handle in which to write CDS sequences
    :param out_genome: file name or handle in which to write genome sequence
    :param identify: function to apply to feature record to get fasta record ID
    :param describe: function to apply to feature record to get fasta record description
    :return:
    """
    logger = logging.getLogger('FastaFromGbk')
    seqs = []
    for record in SeqIO.parse(genbank, 'genbank'):
        seq = SeqRecord(record.seq, id=os.path.splitext(os.path.basename(genbank))[0],
                        description='')
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS':
                    if 'pseudogene' in feature.qualifiers:
                        seq_record = SeqRecord(
                            translate(feature.extract(record.seq), table=11, to_stop=True),
                            id=identify(feature),
                            description=describe(feature))
                    else:
                        seq_record = SeqRecord(
                            Seq(feature.qualifiers['translation'][0]),
                            id=identify(feature),
                            description=describe(feature))
                    seqs.append(seq_record)
    SeqIO.write(seqs, out_cds, 'fasta')
    SeqIO.write(seq, out_genome, 'fasta')
    return

#
# Helper functions for subset_fasta() that can be used for its `match` argument
#
def is_unannotated(name):
    return name.startswith('MTB')

def is_reference(name):
    return not name.startswith(('MTB','L_','L2_'))

def subset_fasta(inseq, outseq, match, identify = lambda _:_):
    """
    write a new fasta file containing only sequences with
    matching record IDs.
    :param infile: str input fasta file name
    :param outseq: str output fasta file name
    :param match: function to apply to record.id for a Boolean result
    :param identify: function to transform record.id prior to matching
    """
    seqs = []
    for record in SeqIO.parse(inseq, 'fasta'):
        record.id = identify(record.id)
        if match(record.id):
            seqs.append(record)
    SeqIO.write(seqs, outseq, 'fasta')

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
