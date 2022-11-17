import os
import subprocess

from Bio import SeqIO


def convert_gbk_to_gff(gbk_filename):
    """
    Runs a subprocess seqret call to convert the input
    gbk_filename to a GFF

    :param gbk_filename: str filename of a Genbank file that needs to be converted to GFF
    :return: None
    """
    seqret_cmd = ['seqret',
                  gbk_filename,
                  '/dev/null',
                  '-oufo',
                  os.path.splitext(gbk_filename)[0] + '.gff',
                  '-feature', '-off', 'gff']
    with open(os.devnull, 'w') as devnull:
        subprocess.call(seqret_cmd, stderr=devnull)


def convert_embl_to_gbk(embl_filename):
    """
    Runs a subprocess seqret call to convert the input
    embl_filename to a GFF

    :param embl_filename: str filename of a EMBL file that needs to be converted to Genbank
    :return: None
    """
    seqret_cmd = ['seqret',
                  embl_filename,
                  os.path.splitext(embl_filename)[0] + '.gbk',
                  '-feature', '-osf', 'genbank']
    with open(os.devnull, 'w') as devnull:
        subprocess.call(seqret_cmd, stderr=devnull)
    return embl_filename.rstrip('embl') + 'gbk'


def convert_gbk_to_embl(genbank_filename):
    """
    Converts a Genbank to EMBL, one EMBL file per Genbank record

    :param genbank_filename: str filename of a EMBL file that needs to be converted to Genbank
    :return: list of embl file names
    """
    embl_files = []
    for record in SeqIO.parse(genbank_filename, "genbank"):
        embl_file = '.'.join([genbank_filename.rstrip('gbk'), record.id, 'embl'])
        SeqIO.write(record, embl_file, "embl")
        embl_files.append(embl_file)
    return embl_files
