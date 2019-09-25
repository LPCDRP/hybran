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
    refdir = args.output + 'temp_references/'
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
    if len(all_embls) >= 30:
        embls = all_embls[0:30]
    else:
        embls = all_embls
    if sorted(embls) != sorted(os.listdir(refdir + '/embls/')):
        logger.info('Getting first ' + str(len(embls)) + ' reference annotations')
        for e in embls:
            gbk = e.split('/')[-1].split('.')[0] + '.gbk'
            gff = e.split('/')[-1].split('.')[0] + '.gff'
            try:
                shutil.copyfile(args.references + e, embl_dir + e)
                shutil.copyfile(args.references + gbk, refdir + gbk)
                shutil.copyfile(args.references + gff, refdir + gff)
            except OSError:
                continue
    return refdir, embl_dir, embls
