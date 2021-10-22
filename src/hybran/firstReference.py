import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_first_reference_proteome(genbank, out_cds, out_genome):
    """

    :param genbank:
    :param out_cds: file handle in which to write CDS sequences
    :param out_genome: file handle in which to write genome sequence
    :return:
    """
    logger = logging.getLogger('ReferenceDatabase')
    logger.info('Creating a reference proteome FASTA for Prokka from ' + genbank)
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
    for s in seqs:
        SeqIO.write(s, out_cds, 'fasta')
    SeqIO.write(seq, out_genome, 'fasta')
    return
