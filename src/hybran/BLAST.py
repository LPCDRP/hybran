from collections import defaultdict
import logging
import os
import subprocess
import tempfile

import multiprocessing
from Bio.SeqRecord import SeqRecord

from functools import partial
from . import config
from .bio import SeqIO


# blast wants to phone home, but we run blast too many times
# and don't want to slow down the program, especially if
# you have no network connection.
os.environ["BLAST_USAGE_REPORT"] = "false"

blast_prog = {
    'n': 'blastn',
    'p': 'blastp',
}

def reference_match(query, subject, seq_ident, seq_covg, identify=lambda _:_, metric='iden', blast_type = "p", strict=False):
    """
    Wrapper to blast(), summarize(), and top_hit() that tells you the bottom line.
    See those functions for descriptions of the main input parameters.
    This function checks for complete hits and, failing that, partial hits.
    :param strict: Boolean whether to strictly enforce all thresholds or allow for low-coverage alignments.
    :returns:
        - result (:py:class:`str`)  -
             name of the top hit (`None` if no suitable match was found)
        - low_covg (:py:class:`bool`) -
             A hit that has either low subject coverage or low query coverage but meets
             remaining thresholds. This indicates whether the query is truncated with
             respect to the reference or vice-versa
        - stats (:py:class:`dict`)  -
             summary dict of blast results from either hits, truncation_signatures, or misses
    """
    
    hits, misses, truncation_signatures = blast(
        query,
        subject,
        seq_ident,
        seq_covg,
        blast_type,
    )
    result = None
    low_covg = False
    hit_dict = None
    if strict:
        principal_lists = [hits]
        scraps = truncation_signatures + misses
    else:
        principal_lists = [hits, truncation_signatures]
        scraps = misses
    # this will be replaced by actual hits later if there are any
    hit_dict = summarize(scraps, identify=identify)

    for hit_list in principal_lists:
        if not hit_list:
            continue
        hit_dict = summarize(hit_list, identify=identify)
        result = top_hit(
            hit_dict,
            metric=metric,
        )
        stats = hit_dict[result]
        scov_pass = stats['scov'] >= seq_covg
        qcov_pass = stats['qcov'] >= seq_covg
        ident_pass = stats['iden'] >= seq_ident
        if ident_pass and any([scov_pass, qcov_pass]) and not all([scov_pass, qcov_pass]):
            low_covg = True
        # if we found a proper hit, don't look for anything else
        break

    return result, low_covg, hit_dict

def top_hit(blast_summary, metric='iden'):
    """
    Get the best hit according to the given metric.

    :param blast_summary: dictionary of blast results as output by summarize()
    :param metric: str inner dictionary key of blast_summary by which to maximize
    :returns: str subject id of the best hit
    """
    return max(blast_summary, key=lambda gene: blast_summary[gene][metric])

def summarize(blast_results, identify = lambda _:_):
    """
    Convert output of blastp() into a dictionary of the following form:
    {id1: {'iden': %-identity, 'scov':scov, 'qcov':qcov}, ...}
    If there are multiple hits against a single subject ID,
    the top-most hit as sorted by blast is retained.

    :param blast_results: list of strings as output by blastp()
    :param identify: function to apply to the subject_id column of
                     the blast results whose value will be used as
                     the outer dictionary key.
    :returns: dictionary of subject IDs to metrics to values that
              defaults to a null hit for missing keys.
    """
    summary = defaultdict(lambda : {
        'iden':0.,
        'qcov':0.,
        'scov':0.,
    })
    for line in reversed(blast_results):
        fields = line.split('\t')
        summary[identify(fields[1])] = {
            'iden':float(fields[2]),
            'qcov':float(fields[14]),
            'scov':float(fields[15]),
        }
    return summary

def blastn(query, subject, seq_ident, seq_covg):
        return blast(query, subject, seq_ident, seq_covg, blast_type = "n")
        
def blastp(query, subject, seq_ident, seq_covg):
        return blast(query, subject, seq_ident, seq_covg, blast_type = "p")

def blast(query, subject, seq_ident, seq_covg, blast_type = "p"):
    """
    Runs BLAST with one sequence as the query and
    the fasta file as the subject.

    :param query: SeqRecord object containing the query sequence/id/description
    :param subject: str subject fasta file name or SeqRecord
    :param seq_ident: sequence identity threshold
    :param seq_covg: alignment coverage threshold
    :param blast_type: nucleotide or protein blast
    :returns:
        - blast_filtered (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits meeting the thresholds
        - blast_rejects (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits falling short of the thresholds
    	- blast_truncation (:py:class:'list') - list of tab-delimited strings corresponding
   					      to BLAST hits where %identity AND either qcov or scov
					      meet the threshold, while the remaining coverage factor fails.

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
    blast_cmd = [
        blast_prog[blast_type],
        '-subject', fa,
        '-outfmt', blast_outfmt,
    ]
    blast_ps = subprocess.run(
        blast_cmd,
        input=str(query.seq),
        text=True,
        check=True,
        capture_output=True,
    )
    stdout = blast_ps.stdout + dummy_hit
    blast_filtered = []
    blast_rejects = []
    blast_truncation = []
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
            elif (identity >= seq_ident) and (qcov >= seq_covg) and (scov <= seq_covg):
                blast_truncation.append('\t'.join(column))
            elif (identity >= seq_ident) and (qcov <= seq_covg) and (scov >= seq_covg):
                blast_truncation.append('\t'.join(column))
            else:
                # Capture the rejected hits
                blast_rejects.append('\t'.join(column))

    return blast_filtered, blast_rejects, blast_truncation


def iterate(fa, seq_list, nproc, seq_ident, seq_covg):
    """
    Runs BLAST function for each query.
    Returns list of all results.
    """
    partial_blast = partial(blastp, subject=fa, seq_ident=seq_ident, seq_covg=seq_covg)
    pool = multiprocessing.Pool(int(nproc))
    hits, misses, truncation_signatures = zip(*pool.map(partial_blast,seq_list))
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
