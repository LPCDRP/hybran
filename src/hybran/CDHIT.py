import logging
from Bio import SeqIO
import subprocess
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_allseq_dict(fa):
    """
    Parses FASTA input, creates dictionary with sequence ID as key and the sequence as the value

    :param fa: str fasta file name
    :returns:
         - fasta_dict (:py:class:`dict`) protein sequence dict by ID
         - fasta_descrip (:py:class:`dict`) sequence descriptions by ID
    """
    fasta_dict={}
    fasta_descrip={}
    for record in SeqIO.parse(fa, "fasta"):
        sequence = str(record.seq)
        fasta_dict[record.id] = sequence
        fasta_descrip[record.id] = record.description
    return fasta_dict, fasta_descrip


def run_cdhit(nproc, input, output, seq_ident, seq_covg):
    """Runs cdhit.

    Coverage flag descriptions from the CD-HIT wiki:
    -aL alignment coverage for the longer sequence
    -aS alignment coverage for the shorter sequence

    CDHIT outputs a fasta file of representative sequences (cdhit_clusters.fasta)
    and a text file of list of clusters (cdhit_clusters.fasta.clstr)."""
    cmd = ['cd-hit',
           '-i', input,
           '-o', output,
           '-c', str(seq_ident * 0.01),
           '-aL',str(seq_covg * 0.01),
           '-aS',str(seq_covg * 0.01),
           '-T', str(nproc),
           '-g', '1',
           '-s', '1',
           '-d', '256']
    out = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return out


def create_reps_dict(in_clusters):
    """
    Creates a dictionary with the representative sequence ID from each cluster as the key
    and the corresponding sequences as the value.

    :param in_clusters: str file name of clusters
    :returns:
         - rep_dict (:py:class:`dict`) cluster members by representatives
         - reps (:py:class:`list`) representatives
    """
    clusters = []
    rep_dict = {}
    rep = ''
    reps = []
    with open(in_clusters, 'r') as f:
        for line in f:
            # At new cluster
            if line.startswith('>'):
                if not rep:
                    continue
                rep_dict[rep] = clusters
                clusters = []
                continue
            # Representatives of cluster
            else:
                if line.rstrip('\n').endswith('*'):
                    rep = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                    reps.append(rep)
                else:
                    gene = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                    clusters.append(gene)
    rep_dict[rep] = clusters
    return rep_dict, reps


def create_reps_fasta(output, reps, rep_seq_id_dict, description_dict):
    """
    Creates a FASTA file that has sequence ID
    and sequence for each representative

    :param output: str output file name
    :param reps: list of representatives
    :param rep_seq_id_dict: dict of representatives and sequences
    :return: None
    """
    seq_list = []
    for id_key in reps:
        record = SeqRecord(Seq(rep_seq_id_dict[id_key]),
                           id=id_key,
                           description=description_dict[id_key])
        seq_list.append(record)
    SeqIO.write(seq_list, output, "fasta")


def run(nproc, fasta, seq_ident,seq_covg, out):
    """
    Executes CDHIT

    :param nproc: int number of processors
    :param fasta: str FASTA file name
    :param seq_ident: float minimum sequence identity
    :param seq_covg: float minimum alignment coverage
    :param out: str output file name
    :return: dict clusters by representative
    """
    logger = logging.getLogger('CDHIT')
    logger.info('Running CDHIT')
    cdhit_ps = run_cdhit(nproc, fasta, out, seq_ident, seq_covg)
    try:
        cdhit_ps.check_returncode()
    except subprocess.CalledProcessError:
        logger.error('CDHIT failed.')
        print(cdhit_ps.stdout, file=sys.stderr)
        exit(1)
    OGdict, description_dict = create_allseq_dict(fasta)
    logger.info('Parsing CDHIT output (' + out + '.clstr)')
    REPdict, rep_list = create_reps_dict(out + '.clstr')
    logger.info('Creating FASTA of CDHIT clusters representatives')
    create_reps_fasta(out, rep_list, OGdict, description_dict)
    return REPdict
