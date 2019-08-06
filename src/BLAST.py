import logging
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing
from functools import partial


def create_raw_seq_list(input_file):
    """
    Parses the fasta file of sequences and
    creates a list containing each header + sequence string.
       """
    seq_list = []
    for record in SeqIO.parse(input_file, 'fasta'):
        seq_list.append([record.id, str(record.seq)])
    return seq_list


def blast(seq_string, fa):
    """
    Runs BLAST with one sequence as the query and
    the fasta file as the subject.
    Returns a list of BLAST runs that
    meet the 95% identity threshold
    """
    blast_outfmt = ' "6 qseqid sseqid pident length mismatch gapopen ' \
                   'qstart qend sstart send evalue bitscore qlen slen"'
    blast_to_all = NcbiblastpCommandline(subject=fa, outfmt=blast_outfmt)
    stdout, stderr=blast_to_all(stdin=seq_string[1])
    blast_filtered = []
    for line in stdout.split('\n'):
        if line:
            column = line.split('\t')
            identity = float(column[2])
            length = float(column[3])
            qlen = float(column[12])
            slen = float(column[13])
            if identity > 95 and length / qlen >= 0.95 and length / slen >= 0.95:
                column[0] = seq_string[0]
                blast_filtered.append('\t'.join(column))
    return blast_filtered


def iterate(fa, seq_list, nproc):
    """
    Runs BLAST function for each query.
    Returns list of all results.
    """
    partial_blast = partial(blast, fa=fa)
    pool = multiprocessing.Pool(nproc)
    list_of_lists = pool.map(partial_blast,seq_list)
    pool.close()
    pool.join()
    all_results_list = map('\n'.join, list_of_lists)
    return all_results_list


def write(all_results_list):
    """
    Joins list of all results into a continuous
    string and writes this string to a new file
    """
    joined_string = '\n'.join(all_results_list)
    with open('blast_results', 'w') as all_v_all:
        all_v_all.write(joined_string)
    return all_v_all


def run_blast(fastafile, nproc):
    logger = logging.getLogger('BLAST')
    logger.info('Running pairwise all-against-all BLAST on ' + fastafile + ' using ' + str(nproc) + ' CPUs')
    seq_string_list = create_raw_seq_list(fastafile)
    all_results_list = iterate(fastafile, seq_string_list, nproc)
    logger.info('Writing BLAST results to blast_results')
    write(all_results_list)
