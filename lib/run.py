import os
import logging
import subprocess

from lib import fastaFromGFF, BLAST, CDHIT, MCL, addEggnogAnnotation, parseClustering


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
        ('ratt' not in os.listdir(isolate) and
         'prokka' not in os.listdir(isolate) and
         'prokka-noreference' not in os.listdir(isolate)) or \
            ('ratt-done' not in os.listdir(isolate + '/ratt/') and
            isolate + '.gbk' not in os.listdir(isolate + '/prokka/') and
            isolate + '.gbk' not in os.listdir(isolate + '/prokka-noreference/')):
        logger.info('Executing RATT and Prokka on ' + isolate)
        try:
            os.mkdir(isolate)
        except OSError:
            pass
        os.chdir(isolate)
        cmd = [script_dir + '/lib/ratt_prokka.sh',
               ref_dir,
               fasta,
               isolate,
               ref_cds,
               cpus]
        subprocess.call(cmd)
        os.chdir(c)


def clustering(annotations, nproc):
    """
    Runs the clustering pipeline which uses CDHIT and MCL
    to cluster orthologous genes. File IO is handled by
    each function

    :param annotations: str directory of all annotations created by annomerge
    :param nproc: str number of processors
    :return: None
    """
    logger = logging.getLogger('ClusterProteins')
    c = os.getcwd()
    try:
        os.mkdir('clustering')
    except OSError:
        pass
    os.chdir('clustering')
    fasta = 'cds_seqs.fasta'
    if 'clustered_proteins' not in os.listdir(os.getcwd()):
        gff_gene_dict = {}
        for d in annotations:
            logger.info('Parsing GFFs in ' + d)
            gff_gene_dict.update(fastaFromGFF.create_fasta(directory=d))
        clusters = CDHIT.cd_hit(nproc=nproc,
                                fasta=fasta,
                                out='cdhit_clusters.fasta')
        BLAST.run_blast(fastafile='cdhit_clusters.fasta',
                        nproc=nproc)
        MCL.run_mcl(in_blast='blast_results',
                    cdhit_clusters=clusters,
                    out_name='clustered_proteins',
                    gene_names=gff_gene_dict)
    os.chdir(c)
    # exit()
    parseClustering.parseClustersUpdateGBKs(gffs=annotations,
                                            clusters='clustering/clustered_proteins')


def eggnog_mapper(script_dir, nproc, emapper_loc):
    """
    Runs the run_emapper.sh which executes emapper.py for
    both diamond and hmmer algorithms. The shell script handles
    all file IO

    :param script_dir: str absolute path to run_emapper.sh
    :param nproc: str number of processors
    :param emapper_loc: str absolute path to eggnog database
    :return: None
    """
    logger = logging.getLogger('OrthologousAnnotation')
    logger.info('Functional annotation with eggnog_mapper')
    try:
        os.mkdir('eggnog-mapper-annotations')
    except OSError:
        pass
    cmd = [script_dir + '/lib/run_emapper.sh',
           nproc,
           emapper_loc]
    subprocess.call(cmd)
    addEggnogAnnotation.update_gbks()
