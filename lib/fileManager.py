import os
import shutil
import logging


def full_path(p):
    """

    :param p:
    :return:
    """
    if not p.starswith('/'):
        return os.getcwd() + p
    return p


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
    if len(all_embls) >= 30:
        embls = all_embls[0:29]
    else:
        embls = all_embls
    if sorted(embls) != sorted(os.listdir(refdir + '/embls/')):
        logger.info('Getting first ' + str(len(embls)) + ' reference annotations')
        for e in embls:
            try:
                shutil.copyfile(args.references + e, embl_dir + e)
                shutil.copyfile(args.references + e, refdir + e.split('.')[0] + '.gbk')
                shutil.copyfile(args.references + e, refdir + e.split('.')[0] + '.gff')
            except OSError:
                continue
    return refdir, embl_dir, embls
