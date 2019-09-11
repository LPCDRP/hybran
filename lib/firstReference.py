import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_first_reference_proteome(genbank):
    """

    :param genbank:
    :return:
    """
    logger = logging.getLogger('ReferenceDatabase')
    logger.info('Creating a reference proteome FASTA for Prokka from ' + genbank)
    seqs = []
    for record in SeqIO.parse(genbank, 'genbank'):
        seq = SeqRecord(record.seq, id=genbank.split('/')[-1].split('.')[0], description='')
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS':
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
    with open('ref_proteome.fasta', 'w') as output:
        for s in seqs:
            SeqIO.write(s, output, 'fasta')
    with open('ref.fasta', 'w') as output:
        SeqIO.write(seq, output, 'fasta')
    return 'ref_proteome.fasta', 'ref.fasta'