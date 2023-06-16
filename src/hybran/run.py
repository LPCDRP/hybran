import os
import logging
import re
import subprocess

from . import extractor, BLAST, CDHIT, MCL, addEggnogAnnotation, parseClustering
from . import config

def ratt_prokka(ref_dir, organism, strain, fasta, ref_cds, gcode, ratt_ttype, prokka_extra_args, script_dir, cpus, qcov):
    """
    Executes RATT and Prokka that resides in ratt_prokka.sh
    File IO is handled by ratt_prokka.sh

    :param ref_dir: str directory that houses the EMBL reference(s)
    :param organism: str genus or binomial name
    :param strain: str strain name
    :param fasta: str FASTA file name that needs to be annotated
    :param ref_cds: str FASTA proteome of the reference for Prokka
    :param gcode: int NCBI genetic code table ID
    :param ratt_ttype: str RATT transfer type
    :param prokka_extra_args: str additional command line flags and arguments for Prokka
    :param script_dir: str absolute path to ratt_prokka.sh
    :param cpus: int number of processors/cpus
    :param qcov: int minimum % query coverage (Prokka doesn't have a way of setting ref coverage)
    :return: None
    """
    logger = logging.getLogger('ProkkaRATTAnnomerge')
    c = os.getcwd()
    isolate = os.path.splitext(os.path.basename(fasta))[0]
    org_components = organism.split()
    if len(org_components) == 2:
        [genus, species] = org_components
        org_flags = f"--genus {genus} --species {species}"
    else:
        org_flags = f"--genus {organism}"
    if strain:
        org_flags += f" --strain {strain}"

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
        cmd = [os.path.join(script_dir,'ratt_prokka.sh'),
               ref_dir,
               fasta,
               isolate,
               ref_cds,
               str(cpus),
               str(qcov),
               str(gcode),
               ratt_ttype,
               f"{org_flags} {prokka_extra_args}",
               ]
        try:
            subprocess.run(
                cmd,
                check=True,
            )
        except subprocess.CalledProcessError:
            logger.error("Could not annotate " + isolate + ".")
            exit(1)
        os.chdir(c)


def clustering(all_genomes, target_genomes, nproc, seq_ident, seq_covg):
    """
    Runs the clustering pipeline which uses CDHIT and MCL
    to cluster orthologous genes. File IO is handled by
    each function

    :param annotations: str directory of all annotations created by annomerge
    :param nproc: int number of processors
    :return: None
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    c = os.getcwd()
    try:
        os.mkdir(c + '/clustering')
    except OSError:
        pass
    try:
        os.mkdir(hybran_tmp_dir + '/clustering')
    except OSError:
        pass
    os.chdir(hybran_tmp_dir + '/clustering')
    fasta = 'cds_seqs.fasta'
    gbk_filenames = [re.sub(r"\.gff$", ".gbk", _) for _ in all_genomes]
    if 'clustered_proteins' not in os.listdir(os.getcwd()):
        gff_gene_dict = {}
        gff_gene_dict.update(extractor.fastaFromGffList(gffs=all_genomes, out_cds=fasta))
        # Run CD-HIT on cds_seqs.fasta
        clusters = CDHIT.run(nproc=nproc,
                             fasta=fasta,
                             seq_ident=seq_ident,
                             seq_covg=seq_covg,
                             out='cdhit_clusters.fasta')
        if 'blast_results' not in os.listdir(hybran_tmp_dir):
            BLAST.run_blast(fastafile='cdhit_clusters.fasta',
                            nproc=nproc,
                            seq_ident=seq_ident,
                            seq_covg=seq_covg)
        if 'clustered_proteins' not in os.listdir(os.getcwd()):
            MCL.run_mcl(in_blast=hybran_tmp_dir + '/blast_results',
                        cdhit_clusters=clusters,
                        out_name='clustered_proteins',
                        gene_names=gff_gene_dict)
    os.chdir(c)
    parseClustering.parseClustersUpdateGBKs(target_gffs=all_genomes,
                                            clusters=hybran_tmp_dir +
                                            '/clustering/clustered_proteins',
                                            genomes_to_annotate=target_genomes,
                                            seq_ident=seq_ident,
                                            seq_covg=seq_covg)


def eggnog_mapper(script_dir, nproc, emapper_loc, ref_tax_ids, ref_gene_dict, temp_dir):
    """
    Runs the run_emapper.sh which executes emapper.py for
    both diamond and hmmer algorithms. The shell script handles
    all file IO

    :param script_dir: str absolute path to run_emapper.sh
    :param nproc: int number of processors
    :param emapper_loc: str absolute path to eggnog database
    :param ref_tax_ids: list of str NCBI taxonomy IDs
    :param ref_gene_dict: dict mapping locus tags to gene names
    :param temp_dir: str path to Hybran temporary directory
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
           temp_dir,
           ','.join(ref_tax_ids),
           ]
    subprocess.call(cmd)
    addEggnogAnnotation.update_gbks(ref_tax_ids, ref_gene_dict)
