import logging
import os
import tempfile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing
from Bio.SeqRecord import SeqRecord

from functools import partial
from . import config


# blast wants to phone home, but we run blast too many times
# and don't want to slow down the program, especially if
# you have no network connection.
os.environ["BLAST_USAGE_REPORT"] = "false"

def summarize(blast_results):
    """
    Convert output of blastp() into a dictionary of the following form:
    {id: [iden,scov,qcov]}
    If there are multiple hits against a single subject ID,
    the top-most hit as sorted by blast is retained.
    """

    return {_[1]: [float(_[2]), float(_[15]), float(_[14])] for _ in
            [line.split('\t') for line in reversed(blast_results)]}

def blastp(query, subject, seq_ident, seq_covg):
    """
    Runs BLAST with one sequence as the query and
    the fasta file as the subject.

    :param query: SeqRecord object containing the query sequence/id/description
    :param subject: str subject fasta file name or SeqRecord
    :param seq_ident: sequence identity threshold
    :param seq_covg: alignment coverage threshold
    :returns:
        - blast_filtered (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits meeting the thresholds
        - blast_rejects (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits falling short of the thresholds
    """
    if isinstance(subject, SeqRecord):
        subject_id = subject.id
        with tempfile.NamedTemporaryFile(
                suffix='.fasta',
                dir=config.hybran_tmp_dir,
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

    # column orders up until bitscore are as expected for mcxdeblast --m9
    # <https://github.com/JohannesBuchner/mcl/blob/5208b974324621f510abb6a29e046a38f6d85f10/src/alien/oxygen/src/mcxdeblast#L259>
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
    stdout += dummy_hit
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

    return blast_filtered, blast_rejects


def iterate(fa, seq_list, nproc, seq_ident, seq_covg):
    """
    Runs BLAST function for each query.
    Returns list of all results.
    """
    partial_blast = partial(blastp, subject=fa, seq_ident=seq_ident, seq_covg=seq_covg)
    pool = multiprocessing.Pool(int(nproc))
    hits, misses = zip(*pool.map(partial_blast,seq_list))
    pool.close()
    pool.join()
    all_results_list = list(map('\n'.join, list(hits)))
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
