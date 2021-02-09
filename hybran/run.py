import os
import logging
import subprocess

from . import fastaFromGFF, BLAST, CDHIT, MCL, addEggnogAnnotation, parseClustering
from . import config

def ratt_prokka(ref_dir, fasta, ref_cds, script_dir, cpus):
    """
    Executes RATT and Prokka that resides in ratt_prokka.sh
    File IO is handled by ratt_prokka.sh

    :param ref_dir: str directory that houses the EMBL reference(s)
    :param fasta: str FASTA file name that needs to be annotated
    :param ref_cds: str FASTA proteome of the reference for Prokka
    :param script_dir: str absolute path to ratt_prokka.sh
    :param cpus: str number of processors/cpus
    :return: None
    """
    logger = logging.getLogger('ProkkaRATTAnnomerge')
    c = os.getcwd()
    isolate = fasta.split('/')[-1].split('.')[0]

    if isolate not in os.listdir(os.getcwd()) or \
        ('ratt' not in os.listdir(isolate) or
         'prokka' not in os.listdir(isolate) or
         'prokka-noreference' not in os.listdir(isolate)) or \
            ('ratt-done' not in os.listdir(isolate + '/ratt/') or
            isolate + '.gbk' not in os.listdir(isolate + '/prokka/') or
            isolate + '.gbk' not in os.listdir(isolate + '/prokka-noreference/')):
        logger.info('Executing RATT and Prokka on ' + isolate)
        try:
            os.mkdir(isolate)
        except OSError:
            pass
        os.chdir(isolate)
        cmd = [os.sep.join([script_dir,'ratt_prokka.sh']),
               ref_dir,
               fasta,
               isolate,
               ref_cds,
               str(cpus)]
        subprocess.call(cmd)
        os.chdir(c)


def clustering(all_genomes, target_genomes, nproc, seq_ident):
    """
    Runs the clustering pipeline which uses CDHIT and MCL
    to cluster orthologous genes. File IO is handled by
    each function

    :param annotations: str directory of all annotations created by annomerge
    :param nproc: str number of processors
    :return: None
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    c = os.getcwd()
    try:
        os.mkdir(hybran_tmp_dir + '/clustering')
    except OSError:
        pass
    os.chdir(hybran_tmp_dir + '/clustering')
    fasta = 'cds_seqs.fasta'
    if 'clustered_proteins' not in os.listdir(os.getcwd()):
        gff_gene_dict = {}
        gff_gene_dict.update(fastaFromGFF.create_fasta(directory=all_genomes))
        clusters = CDHIT.run(nproc=nproc,
                             fasta=fasta,
                             seq_ident=seq_ident,
                             out='cdhit_clusters.fasta')
        if 'blast_results' not in os.listdir(hybran_tmp_dir):
            BLAST.run_blast(fastafile='cdhit_clusters.fasta',
                            nproc=nproc)
        if 'clustered_proteins' not in os.listdir(os.getcwd()):
            MCL.run_mcl(in_blast=hybran_tmp_dir + '/blast_results',
                        cdhit_clusters=clusters,
                        out_name='clustered_proteins',
                        gene_names=gff_gene_dict)
    os.chdir(c)
    parseClustering.parseClustersUpdateGBKs(target_gffs=all_genomes,
                                            clusters=hybran_tmp_dir +
                                            '/clustering/clustered_proteins',
                                            genomes_to_annotate=target_genomes)


def eggnog_mapper(script_dir, nproc, emapper_loc, temp_dir):
    """
    Runs the run_emapper.sh which executes emapper.py for
    both diamond and hmmer algorithms. The shell script handles
    all file IO

    :param script_dir: str absolute path to run_emapper.sh
    :param nproc: str number of processors
    :param emapper_loc: str absolute path to eggnog database
    :param temp_dir: str path to hybran temporary directory
    :return: None
    """

    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('OrthologousAnnotation')
    logger.info('Functional annotation with eggnog_mapper')
    try:
        os.mkdir(hybran_tmp_dir + '/eggnog-mapper-annotations')
    except OSError:
        pass
    cmd = [os.sep.join([script_dir, 'run_emapper.sh']),
           str(nproc),
           emapper_loc,
           temp_dir]
    subprocess.call(cmd)
    addEggnogAnnotation.update_gbks(script_dir)
