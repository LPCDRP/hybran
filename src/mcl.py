import subprocess
import logging
import networkx as nx


def execute_mcl(blast):
    mcxdeblast_cmd = ['mcxdeblast',
                      '-m9',
                      '--score=r',
                      '--line-mode=abc',
                      '--out=mcxdeblast_results',
                      blast]
    mcl_cmd = ['mcl',
               'mcxdeblast_results',
               '--abc',
               '-I', '1.5',
               '-o', 'mcl',
               '-q', 'x',
               '-V', 'all']
    mcxdeblast_out = subprocess.Popen(mcxdeblast_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    mcl_out = subprocess.Popen(mcl_cmd,
                               stdout=subprocess.PIPE,)


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
    execute_mcl(blast=in_blast)
    logger.info('Re-inflating CDHIT clusters in MCL clusters')
    output = inflate_mcl_clusters(mcl_output='mcl',
                                  cdhit_groups=cdhit_clusters,
                                  gene_key=gene_names)
    logger.info('Writing final output ' + out_name)
    writer(output, out_name)
