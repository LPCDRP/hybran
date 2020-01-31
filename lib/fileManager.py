import os
import shutil
import logging
import converter


def full_path(p):
    """

    :param p:
    :return:
    """
    if not p.startswith('/'):
        return os.getcwd() + '/' + p
    return p


def remove_file(f_list):
    for f in f_list:
        try:
            os.remove(f)
        except OSError:
            pass


def remove_dir(d_list):
    for d in d_list:
        try:
            shutil.rmtree(d)
        except OSError:
            pass


def ratt_references(args):
    logger = logging.getLogger('RATTReferences')
    refdir = args.references + 'temp_references/'
    embl_dir = refdir + '/embls/'

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
            if embl_count <= 30:
                shutil.copyfile(args.references + e, embl_dir + e)
                embls.append(e)
            shutil.copyfile(args.references + gbk, refdir + gbk)
            shutil.copyfile(args.references + gff, refdir + gff)
            embl_count += 1
        except OSError:
            continue
    return refdir, embl_dir, embls
