import logging
import os
import subprocess
import sys
import time

from . import config


def execute_mcxdeblast(blast):
    hybran_tmp_dir = config.hybran_tmp_dir
    mcxdeblast_cmd = ['mcxdeblast',
                      '-m9',
                      '--score=r',
                      '--line-mode=abc',
                      '--out=' + hybran_tmp_dir + '/mcxdeblast_results',
                      blast]
    mcxdeblast_out = subprocess.run(
        mcxdeblast_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return mcxdeblast_out


def inflate_mcl_clusters(mcl_output, cdhit_groups, gene_key):
    """Produces the contents of the clustered_proteins file, which is
    some combination of the MCL output and the CD-HIT ouput."""
    inflated_output = []
    with open(mcl_output, 'r') as mcl:
        for line in mcl:
            names = line.rstrip('\n').split('\t')
            try:
                inflated_cluster = names + list(set(sorted([m for n in names for m in cdhit_groups[n]])))
            except KeyError:
                inflated_cluster = names
            if gene_key:
                gene_name = gene_key[inflated_cluster[0]]
            else:
                gene_name = inflated_cluster[0]
            inflated_output.append([gene_name + ': ' + inflated_cluster[0]] + inflated_cluster[1:])
    return inflated_output


def writer(lines, output_name):
    with open(output_name, 'w') as out:
        for line in lines:
            out.write('\t'.join(line) + '\n')


def run_mcl(in_blast, cdhit_clusters, out_name, gene_names, inflation):
    hybran_tmp_dir = config.hybran_tmp_dir    
    logger = logging.getLogger('MCL')
    logger.info('Running MCL')
    mcx_ps = execute_mcxdeblast(blast=in_blast)
    try:
         mcx_ps.check_returncode()
         logger.debug('\n' + mcx_ps.stdout)
    except subprocess.CalledProcessError:
        logger.error('\n' + mcx_ps.stdout)
        logger.error('mcxdeblast failed.')
    time.sleep(5)
    mcl_cmd = [
        'mcl',
        os.path.join(hybran_tmp_dir, 'mcxdeblast_results'),
        '--abc',
        '-I', str(inflation),
        '-o', os.path.join(hybran_tmp_dir, 'mcl'),
        '-q', 'x',
        '-V', 'all',
    ]
    mcl_ps = subprocess.run(
        mcl_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        mcl_ps.check_returncode()
    except subprocess.CalledProcessError:
        logger.error('\n' + mcl_ps.stdout)
        logger.error('mcl failed')
    time.sleep(5)
    output = inflate_mcl_clusters(mcl_output=hybran_tmp_dir + '/mcl',
                                  cdhit_groups=cdhit_clusters,
                                  gene_key=gene_names)
    logger.info('Writing final output ' + out_name)
    # Write clustered_proteins file
    writer(output, out_name)
