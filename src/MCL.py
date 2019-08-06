import subprocess
import logging
import time


def execute_mcxdeblast(blast):

    mcxdeblast_cmd = ['mcxdeblast',
                      '-m9',
                      '--score=r',
                      '--line-mode=abc',
                      '--out=mcxdeblast_results',
                      blast]
    mcxdeblast_out = subprocess.Popen(mcxdeblast_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)


def execute_mcl():
    mcl_cmd = ['mcl',
               'mcxdeblast_results',
               '--abc',
               '-I', '1.5',
               '-o', 'mcl',
               '-q', 'x',
               '-V', 'all']
    mcl_out = subprocess.Popen(mcl_cmd,
                               stdout=subprocess.PIPE)
    return mcl_out


def inflate_mcl_clusters(mcl_output, cdhit_groups, gene_key):
    inflated_output = []
    with open(mcl_output, 'r') as mcl:
        for line in mcl:
            names = line.rstrip('\n').split('\t')
            try:
                inflated_cluster = sorted([cdhit_groups[n] for n in names][0] + names)
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


def run_mcl(in_blast, cdhit_clusters, out_name, gene_names):
    logger = logging.getLogger('MCL')
    logger.info('Running MCL')
    execute_mcxdeblast(blast=in_blast)
    time.sleep(5)
    mcl_stdout = execute_mcl()
    time.sleep(5)
    logger.info('Re-inflating CDHIT clusters in MCL clusters')
    output = inflate_mcl_clusters(mcl_output='mcl',
                                  cdhit_groups=cdhit_clusters,
                                  gene_key=gene_names)
    logger.info('Writing final output ' + out_name)
    writer(output, out_name)
