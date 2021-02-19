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


def ratt_references(args):
    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('RATTReferences')
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
    all_embls = [i for i in os.listdir(args.references) if i.endswith('.embl')]
    if not all_embls:
        for i in os.listdir(args.references):
            if i.endswith('.gbk'):
                converter.convert_gbk_to_embl(args.references + i)
    embls = []
    embl_count = 0
    logger.info('Getting all reference annotations')
    for e in all_embls:
        gbk = e.split('/')[-1].split('.')[0] + '.gbk'
        gff = e.split('/')[-1].split('.')[0] + '.gff'
        try:
            shutil.copyfile(args.references + e, embl_dir + e)
            embls.append(e)
            shutil.copyfile(args.references + gbk, refdir + gbk)
            shutil.copyfile(args.references + gff, refdir + gff)
            embl_count += 1
        except OSError:
            continue
    return refdir, embl_dir, embls
