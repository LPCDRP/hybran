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

presets = {
    'thorough': [
        # comp_based_stats described in:
        # Yu & Altschul, 2005 doi:10.1093/bioinformatics/bti070
        # Reason for this setting:
        # https://github.com/ncbi/blast_plus_docs/issues/23
        '-comp_based_stats', '0',
        '-seg', 'no',
    ],
    'conservative': [
        '-comp_based_stats', '2',
        # Moreno-Hagelsieb & Latimer, 2008 doi:10.1093/bioinformatics/btm585
        #'-use_sw_tback', # authors said this doesn't add much improvement but slows things down signicantly
        '-seg', 'yes',
        '-soft_masking', 'true',
    ],
}


def reference_match(
        query,
        subject,
        cutoff=config.cnf.blast.min_identity,
        seq_covg=config.cnf.blast.min_coverage,
        identify=lambda _:_,
        metric='iden',
        blast_type="p",
        strict=False,
        preset='thorough',
):
    """
    Wrapper to blast(), summarize(), and top_hit() that tells you the bottom line.
    See those functions for descriptions of the main input parameters.
    This function checks for complete hits and, failing that, partial hits.
    :param query: SeqRecord a single gene to search
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
        query=query,
        subject=subject,
        metric=metric,
        cutoff=cutoff,
        seq_covg=seq_covg,
        blast_type=blast_type,
        blast_extra_args=presets[preset],
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
    hit_dict = summarize(scraps, identify=identify)[query.id]

    for hit_list in principal_lists:
        if not hit_list:
            continue
        hit_dict = summarize(hit_list, identify=identify)[query.id]
        result = top_hit(
            hit_dict,
            metric=metric,
        )
        stats = hit_dict[result]
        scov_pass = stats['scov'] >= seq_covg
        qcov_pass = stats['qcov'] >= seq_covg
        low_covg = (
            stats[metric] >= cutoff
            and any([scov_pass, qcov_pass])
            and not all([scov_pass, qcov_pass])
        )
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

def top_hits(blast_summary, metric='iden'):
    """
    Get the best hits according to the given metric, those with scores
    >= 99% of the highest absolute score.
    See Wolf & Koonin, 2012 (doi:10.1093/gbe/evs100)
    """
    top_score = max([blast_summary[gene][metric] for gene in blast_summary])
    return sorted(
        [gene for gene in blast_summary if blast_summary[gene][metric] >= 0.99 * top_score],
        key=lambda gene: blast_summary[gene][metric],
        reverse=True,
    )

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
    summary = defaultdict(lambda: defaultdict(lambda : {
        'iden':0.,
        #'evalue': 1e10,
        'bitscore':0.,
        'qcov':0.,
        'scov':0.,
    }))
    for line in reversed(blast_results):
        fields = line.split('\t')
        summary[identify(fields[0])][identify(fields[1])] = {
            'iden':float(fields[2]),
            #'evalue': float(fields[10]),
            'bitscore': float(fields[11]),
            'qcov':float(fields[14]),
            'scov':float(fields[15]),
        }

    return summary

def blastn(query, subject, metric, cutoff, seq_covg, nproc=1, blast_extra_args=None, screen=True):
        return blast(query, subject, metric, cutoff, seq_covg, blast_type="n", nproc=nproc, blast_extra_args=blast_extra_args, screen=screen)
        
def blastp(query, subject, metric, cutoff, seq_covg, nproc=1, blast_extra_args=None, screen=True):
        return blast(query, subject, metric, cutoff, seq_covg, blast_type="p", nproc=nproc, blast_extra_args=blast_extra_args, screen=screen)

def blast(
        query,
        subject,
        metric,
        cutoff,
        seq_covg,
        blast_type="p",
        nproc=1,
        blast_extra_args=None,
        screen=True,
):
    """
    Runs BLAST

    :param query: str query fasta file name or SeqRecord object containing the query sequence/id/description
    :param subject: str subject fasta file name or SeqRecord
    :param metric: str name of parameter to use as primary threshold (in addition to alignment coverage)
    :param cutoff: float threshold for selected metric
    :param seq_covg: alignment coverage threshold
    :param blast_type: nucleotide or protein blast
    :param nproc: int number of threads to use for blast
    :param blast_extra_args: list of str tokens to pass to blast executable
    :param screen: bool whether to apply thresholds/filtering to results
    :returns:
        - blast_filtered (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits meeting the thresholds
        - blast_rejects (:py:class:`list`) - list of tab-delimited strings corresponding
                                              to BLAST hits falling short of the thresholds
    	- blast_truncation (:py:class:'list') - list of tab-delimited strings corresponding
   					      to BLAST hits where %identity AND either qcov or scov
					      meet the threshold, while the remaining coverage factor fails.

    """
    # prevent a mutable variable (list) being defined in the function header
    if blast_extra_args is None:
        blast_extra_args = []
    # prevent changes to the argument (which is usually the preset list) from persisting
    else:
        blast_extra_args = blast_extra_args.copy()

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

    if not isinstance(query, SeqRecord):
        blast_extra_args += ['-query', query]
        blast_stdin = None
        count = 1
        for rec in SeqIO.parse(query,'fasta'):
            if count > 1:
                query_id = ''
                break
            query_id = rec.id
            count += 1
    else:
        query_id = query.id
        blast_stdin = str(query.seq)

    # column orders up until bitscore are as expected for mcxdeblast --m9
    # <https://github.com/JohannesBuchner/mcl/blob/5208b974324621f510abb6a29e046a38f6d85f10/src/alien/oxygen/src/mcxdeblast#L259>
    blast_outfmt = "6 qseqid sseqid pident length mismatch gapopen " \
                   "qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    # for a single alignment, add a dummy hit with 0% identity so that 0-thresholds are valid
    # when there are no hits at all.
    if subject_id:
        dummy_hit = '\t'.join([query_id, subject_id,'0','0','0','0',
                               '','','','','1','0','1','1','',''])+'\n'
    else:
        dummy_hit = ''
    blast_cmd = [
        blast_prog[blast_type],
        '-subject', fa,
        '-outfmt', blast_outfmt,
        '-num_threads', str(nproc),
    ] + blast_extra_args
    blast_ps = subprocess.run(
        blast_cmd,
        input=blast_stdin,
        text=True,
        check=True,
        capture_output=True,
    )

    if not screen:
        return blast_ps.stdout.split('\n'), [], []

    stdout = blast_ps.stdout + dummy_hit
    blast_filtered = []
    blast_rejects = []
    blast_truncation = []
    for line in stdout.split('\n'):
        if not line:
            continue
        values = {}

        column = line.split('\t')
        values['iden'] = float(column[2])
        length = float(column[3])
        values['bitscore'] = float(column[11])
        qlen = float(column[12])
        slen = float(column[13])
        qseq, sseq = column[14:16]
        qcov = (len(qseq.replace("-","")) / qlen) * 100
        scov = (len(sseq.replace("-","")) / slen) * 100
        column[14:16] = str(qcov), str(scov)
        if blast_stdin is not None:
            column[0] = query.id
        if (values[metric] >= cutoff) and (qcov >= seq_covg) and (scov >= seq_covg):
            blast_filtered.append('\t'.join(column))
        elif (values[metric] >= cutoff) and (qcov >= seq_covg) and (scov <= seq_covg):
            blast_truncation.append('\t'.join(column))
        elif (values[metric] >= cutoff) and (qcov <= seq_covg) and (scov >= seq_covg):
            blast_truncation.append('\t'.join(column))
        else:
            # Capture the rejected hits
            blast_rejects.append('\t'.join(column))

    return blast_filtered, blast_rejects, blast_truncation


def all_vs_all(
        fastafile,
        nproc,
        metric='iden',
        cutoff=config.cnf.blast.min_identity,
        seq_covg=config.cnf.blast.min_coverage,
        blast_preset='thorough',
):
    logger = logging.getLogger('BLAST')
    logger.info('Running pairwise all-against-all BLAST on ' + fastafile + ' using ' + str(nproc) + ' CPUs')

    all_results, _, _ = blastp(
        query=fastafile,
        subject=fastafile,
        metric=metric,
        cutoff=cutoff,
        seq_covg=seq_covg,
        nproc=nproc,
        blast_extra_args=presets[blast_preset],
        screen=True,
    )

    logger.info('Writing BLAST results to blast_results in Hybran temporary '
                'directory.')
    hybran_tmp_dir = config.hybran_tmp_dir
    with open(hybran_tmp_dir + '/blast_results', 'w') as all_v_all:
        all_v_all.write('\n'.join(all_results))

def bidirectional_best_hit(
        query,
        subject,
        metric='iden',
        cutoff=config.cnf.blast.min_identity,
        min_seq_covg=config.cnf.blast.min_coverage,
        identify=lambda _:_,
        nproc=1,
        blast_type="p",
):
    """
    Identify bidirectional BLAST best hits between a group of query sequences and a group of reference sequences.

    :param query: str fasta file name containing query sequences
    :param subject: str fasta file name containing reference sequences
    :param metric: str name of parameter to use as the primary threshold (in addition to alignment coverage)
    :param cutoff: float threshold for the selected metric
    :param min_seq_covg: float minimum percent sequence alignment coverage
    :param identify:
      function to apply to extract the desired sequence names from the blast results.
      passed to summarize()
    """
    bbh_results = defaultdict(lambda :None)
    blast_extra_args = presets['conservative']

    qry2ref_hits, _, _ = blast(
        query,
        subject,
        metric=metric,
        cutoff=cutoff,
        seq_covg=min_seq_covg,
        blast_type=blast_type,
        nproc=nproc,
        blast_extra_args=blast_extra_args,
    )
    qry2ref_hit_dict = summarize(qry2ref_hits, identify=identify)
    ref2qry_hits, _, _ = blast(
        subject,
        query,
        metric=metric,
        cutoff=cutoff,
        seq_covg=min_seq_covg,
        blast_type=blast_type,
        nproc=nproc,
        blast_extra_args=blast_extra_args,
    )
    ref2qry_hit_dict = summarize(ref2qry_hits, identify=identify)
    ref2qry_top_hits = defaultdict(lambda : [])
    ref2qry_top_hits.update({
        ref_gene:top_hits(ref_hit_summary)
        for ref_gene, ref_hit_summary in ref2qry_hit_dict.items()
    })

    for query_gene in qry2ref_hit_dict:
        bbh_match = None
        score_to_beat = 0.
        qry_top_hits = top_hits(qry2ref_hit_dict[query_gene])
        for ref_gene in qry_top_hits:
            if (
                    query_gene in ref2qry_top_hits[ref_gene]
                    and qry2ref_hit_dict[query_gene][ref_gene][metric] > score_to_beat
            ):
                bbh_match = ref_gene
                score_to_beat = qry2ref_hit_dict[query_gene][ref_gene][metric]
        bbh_results[query_gene] = bbh_match

    return bbh_results, qry2ref_hit_dict
