import glob
import os
import shutil
import logging

from . import converter
from . import config
from . import designator
from .bio import SeqIO


exts = dict(
    fasta=['fna', 'fasta', 'fa'],
    genbank=['gbk','gbf','gb', 'gbff'],
)

def file_list(glist, file_type="fasta"):
    """
    parse a genome list argument, which may contain individual file names, a directory,
    or a file of file names, and return a list of absolute paths to fasta files
    :param glist: a list of strings
    :param file_type: str file format to look for. fasta or genbank
    :return genomes: a list of absolute paths
    """
    if os.path.isdir(glist[0]):
        genomes = []
        for ext in exts[file_type]:
            genomes += glob.glob(os.path.join(glist[0],'*.'+ext))
        genomes = [os.path.abspath(_) for _ in genomes]
    elif len(glist) > 1:
        genomes = [os.path.abspath(g) for g in glist]
    elif any([glist[0].endswith('.'+ext) for ext in exts[file_type]]):
        genomes = [os.path.abspath(glist[0])]
    else:
        genomes = []
        with open(glist[0], 'r') as genome:
            for line in genome:
                genomes.append(line.rstrip('\n'))
    return genomes

def prepare_references(references):
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
    for i in references:
        # we will be using the references exclusively from our temporary directory
        # genbank
        ref_file = os.path.basename(i)
        ref_id = os.path.splitext(ref_file)[0]
        gbk = os.path.join(refdir, ref_file)
        revised_records = []
        for record in SeqIO.parse(i, "genbank"):
            ref_contig_id = '.'.join([ref_id, record.id])
            for f in record.features:
                if ( 'locus_tag' in f.qualifiers.keys()
                     and 'gene' not in f.qualifiers.keys()
                ):
                    f.qualifiers['gene'] = [ f.qualifiers['locus_tag'][0] ]
                # Since RATT doesn't tell us which annotation came from which reference,
                # we'll mark our private copies of the input references so we can trace
                # them after the transfer.
                designator.append_qualifier(
                    f.qualifiers,
                    'note',
                    f"HYBRANSOURCE:{ref_contig_id}",
                )
            revised_records.append(record)
        SeqIO.write(revised_records, gbk, "genbank")

        # embl - RATT needs them in their own directory, but they're sometimes looked for in
        #        the top-level folder.
        embl_recs = converter.convert_gbk_to_embl(gbk)
        for embl in embl_recs:
            os.symlink(embl, os.path.join(embl_dir, os.path.basename(embl)))
        # gff
        converter.convert_gbk_to_gff(gbk)

    embls = glob.glob(os.path.join(refdir, '*.embl'))

    return refdir, embl_dir, embls
