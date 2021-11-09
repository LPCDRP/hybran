import logging
import os
import tempfile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing
from Bio.SeqRecord import SeqRecord

from functools import partial
from . import config


def summarize(blast_results):
    """
    Convert output of blastp() into a dictionary of the following form:
    {id: [iden,scov,qcov]}
    """

    # TODO - this uses the same logic as what was in annomerge, but
    #        it has the problem that multiple hits for the same subject
    #        will get overwritten by the last one in the blast output list.
    return {_[1]: [float(_[2]), float(_[15]), float(_[14])] for _ in
            [line.split('\t') for line in blast_results]}

def blastp(query, subject, seq_ident, seq_covg):
    """
    Runs BLAST with one sequence as the query and
    the fasta file as the subject.

    :param query: SeqRecord object containing the query sequence/id/description
    :param subject: str subject fasta file name or SeqRecord
    :param seq_ident: sequence identity threshold
    :param seq_covg: alignment coverage threshold
    :return: list of tab-delimited strings corresponding to BLAST hits meeting the thresholds
    """
    if isinstance(subject, SeqRecord):
        subject_id = subject.id
        with tempfile.SpooledTemporaryFile(
                suffix='.fasta',
                dir=os.path.join(config.hybran_tmp_dir,'seqs'),
                delete=False,
                mode='w',
        ) as fa_handle:
            fa = fa_handle.name
            SeqIO.write(subject, fa, "fasta")
    else:
        fa = subject
        # quickly check how many sequences are in here. if it's just one,
        # we want to pad the output.
        count = 1
        for rec in SeqIO.parse(fa,'fasta'):
            if count > 1:
                subject_id = ''
                break
            subject_id = rec.id
            count += 1

    blast_outfmt = "6 qseqid sseqid pident length mismatch gapopen " \
                   "qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    # for a single alignment, add a dummy hit with 0% identity so that 0-thresholds are valid
    # when there are no hits at all.
    if subject_id:
        dummy_hit = '\t'.join([query.id, subject_id,'0','0','0','0',
                               '','','','','1','0','1','1','',''])+'\n'
    else:
        dummy_hit = ''
    blast_to_all = NcbiblastpCommandline(subject=fa, outfmt=blast_outfmt)
    stdout, stderr=blast_to_all(stdin=str(query.seq))
    stdout = dummy_hit + stdout
    blast_filtered = []
    blast_rejects = []
    for line in stdout.split('\n'):
        if line:
            column = line.split('\t')
            identity = float(column[2])
            length = float(column[3])
            qlen = float(column[12])
            slen = float(column[13])
            qseq, sseq = column[14:16]
            qcov = (len(qseq.replace("-","")) / qlen) * 100
            scov = (len(sseq.replace("-","")) / slen) * 100
            column[14:16] = str(qcov), str(scov)
            column[0] = query.id
            if (identity >= seq_ident) and (qcov >= seq_covg) and (scov >= seq_covg):
                blast_filtered.append('\t'.join(column))
            else:
                # Capture the rejected hits
                blast_rejects.append('\t'.join(column))

    # Append the rejected hits
    with open('../blast_rejects', 'a') as f:
        for line in blast_rejects:
            f.write(str(line) + '\n')

    return blast_filtered


def iterate(fa, seq_list, nproc, seq_ident, seq_covg):
    """
    Runs BLAST function for each query.
    Returns list of all results.
    """
    partial_blast = partial(blastp, subject=fa, seq_ident=seq_ident, seq_covg=seq_covg)
    pool = multiprocessing.Pool(int(nproc))
    list_of_lists = pool.map(partial_blast,seq_list)
    pool.close()
    pool.join()
    all_results_list = list(map('\n'.join, list_of_lists))
    return all_results_list


def write(all_results_list):
    """
    Joins list of all results into a continuous
    string and writes this string to a new file
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    joined_string = '\n'.join(all_results_list)
    with open(hybran_tmp_dir + '/blast_results', 'w') as all_v_all:
        all_v_all.write(joined_string)
    return all_v_all


def run_blast(fastafile, nproc, seq_ident, seq_covg):
    logger = logging.getLogger('BLAST')
    logger.info('Running pairwise all-against-all BLAST on ' + fastafile + ' using ' + str(nproc) + ' CPUs')
    seq_string_list = SeqIO.parse(fastafile, 'fasta')
    all_results_list = iterate(fastafile, seq_string_list, nproc, seq_ident=seq_ident, seq_covg=seq_covg)
    logger.info('Writing BLAST results to blast_results in Hybran temporary '
                'directory.')
    write(all_results_list)
