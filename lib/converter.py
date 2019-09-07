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
                  gbk_filename.split('.')[0] + '.gff',
                  '-feature', '-osf', 'gff']
    subprocess.Popen(seqret_cmd, stderr=os.devnull)
