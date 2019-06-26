
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
    

def run_cdhit(nproc, input, output):
    """Runs cdhit with a sequence identity threshold of 0.98."""
    cmd = ['cdhit',
           '-i', input,
           '-o', output,
           '-c', '0.98',
           '-T', str(nproc),
           '-g', '1',
           '-s', '1',
           '-d', '256']
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE)


def find_reps(clstr_file):
    """Parses the clstr file (splits on white space) that is now in the current directory 
    after running cdhit. If the last column is a "*" then it is the representative sequence.
    Loop over each line in the file and put the sequence ID's for the representative sequences into a list."""
    rep_seq_id_list = []
    with open(clstr_file, 'r') as clstrs:
        for line in clstrs:
            info = line.strip().split()
            indicator = info[-1]
            if indicator == '*':
                id1 = info[2]
                partial_id = id1[1:-3]
                rep_seq_id_list.append(partial_id)
    return rep_seq_id_list


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
    return rep_dict, reps


def create_reps_fasta(output, reps, rep_seq_id_dict):
    """Creates a FASTA file that has sequence ID and sequence for each representative"""
    seq_list = []
    for id_key in reps:
        record = SeqRecord(Seq(rep_seq_id_dict[id_key]), id=id_key, description="")
        seq_list.append(record)
    SeqIO.write(seq_list, output, "fasta")



def cd_hit(nproc, fasta, out):
    run_cdhit(nproc, fasta, out)
    OGdict = create_allseq_dict(fasta)
    REPdict, rep_list = create_reps_dict(out + ".clstr")
    create_reps_fasta(out, rep_list, OGdict)
    return REPdict
