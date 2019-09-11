import os
import subprocess


def convert_gbk_to_gff(gbk_filename):
    """
    Runs a subprocess seqret call to convert the input
    gbk_filename to a GFF

    :param gbk_filename: str filename of a Genbank file that needs to be converted to GFF
    :return: None
    """
    seqret_cmd = ['seqret',
                  gbk_filename,
                  gbk_filename.split('.')[-2] + '.gff',
                  '-feature', '-osf', 'gff']
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
                  embl_filename.rstrip('embl') + 'gbk',
                  '-feature', '-osf', 'genbank']
    with open(os.devnull, 'w') as devnull:
        subprocess.call(seqret_cmd, stderr=devnull)
    return embl_filename.rstrip('embl') + 'gbk'
