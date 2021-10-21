import glob
import os
import shutil
import logging
from . import converter
from . import config

def full_path(p):
    """

    :param p:
    :return:
    """
    if not p.startswith('/'):
        return os.getcwd() + '/' + p
    return p

def genome_list(glist):
    """
    parse a genome list argument, which may contain individual file names, a directory,
    or a file of file names, and return a list of absolute paths to fasta files
    :param glist: a list of strings
    :return genomes: a list of absolute paths
    """
    fasta_exts = ['fna', 'fasta', 'fa']
    if os.path.isdir(glist[0]):
        genomes = []
        for ext in fasta_exts:
            genomes += glob.glob(os.path.join(glist[0],'*.'+ext))
        genomes = [os.path.abspath(_) for _ in genomes]
    elif len(glist) > 1:
        genomes = [os.path.abspath(g) for g in glist]
    elif any([glist[0].endswith('.'+ext) for ext in fasta_exts]):
        genomes = [os.path.abspath(glist[0])]
    else:
        genomes = []
        with open(glist[0], 'r') as genome:
            for line in genome:
                genomes.append(line.rstrip('\n'))
    return genomes

def prepare_references(args):
    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('PrepareReferences')
    refdir = hybran_tmp_dir + '/temp_references/'
    embl_dir = refdir + 'embls/'

    try:
        os.mkdir(refdir)
    except OSError:
        pass
    try:
        os.mkdir(embl_dir)
    except OSError:
        pass
    embls = []
    embl_count = 0
    logger.info('Getting all reference annotations')
    for i in glob.iglob(os.path.join(args.references, '*.gbk')):
        # we will be using the references exclusively from our temporary directory
        # genbank
        gbk = os.path.join(refdir, os.path.basename(i))
        os.symlink(i, gbk)
        # embl - RATT needs them in their own directory, but they're sometimes looked for in
        #        the top-level folder.
        embl = converter.convert_gbk_to_embl(gbk)
        os.symlink(embl, os.path.join(embl_dir, os.path.basename(embl)))
        # gff
        converter.convert_gbk_to_gff(gbk)

    embls = glob.glob(os.path.join(refdir, '*.embl'))

    return refdir, embl_dir, embls
