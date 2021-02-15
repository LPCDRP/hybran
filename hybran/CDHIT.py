import logging
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_allseq_dict(fa):
    """Parses FASTA input, creates dictionary with sequence ID as key and the sequence as the value"""
    fasta_dict={}
    for record in SeqIO.parse(fa, "fasta"):
        sequence = str(record.seq)
        fasta_dict[record.id] = sequence
    return fasta_dict
    

def run_cdhit(nproc, input, seq_ident, output):
    """Runs cdhit with a sequence identity threshold of 0.95."""
    cmd = ['cd-hit',
           '-i', input,
           '-o', output,
           '-c', str(seq_ident),
           '-T', str(nproc),
           '-g', '1',
           '-s', '1',
           '-d', '256']
    out = subprocess.run(cmd, stdout=subprocess.PIPE)
    return out


def create_reps_dict(in_clusters):
    """Creates a dictionary with the representative sequence ID from each cluster as the key
    and the corresponding sequences as the value."""
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


def create_reps_fasta(output, reps, rep_seq_id_dict):
    """Creates a FASTA file that has sequence ID and sequence for each representative"""
    seq_list = []
    for id_key in reps:
        record = SeqRecord(Seq(rep_seq_id_dict[id_key]), id=id_key, description="")
        seq_list.append(record)
    SeqIO.write(seq_list, output, "fasta")


def run(nproc, fasta, seq_ident, out):
    logger = logging.getLogger('CDHIT')
    logger.info('Running CDHIT')
    run_cdhit(nproc, fasta, seq_ident, out)
    OGdict = create_allseq_dict(fasta)
    logger.info('Parsing CDHIT output (' + out + '.clstr)')
    REPdict, rep_list = create_reps_dict(out + '.clstr')
    logger.info('Creating FASTA of CDHIT clusters representatives')
    create_reps_fasta(out, rep_list, OGdict)
    return REPdict
