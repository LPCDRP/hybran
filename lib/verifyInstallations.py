import os
import logging


def verify_installations(eggnog_databases):
    """

    :param eggnog_databases:
    :return:
    """
    logger = logging.getLogger('VerifyInstallations')

    # Verify the RATT_HOME variable is set
    if 'RATT_HOME' not in os.environ.keys():
        logger.error('RATT_HOME variable is not defined. '
                     'Please export RATT_HOME=[location of downloaded ratt executable/config')
        exit(EnvironmentError)
    # Verify eggnog-mapper/ is present in given directory
    if not os.path.isdir(eggnog_databases):
        logger.error('There are no databases in ' + eggnog_databases + '. Location of eggnog databases is required.'
                                                                       '\nIf the databases have not been downloaded, '
                                                                       'execute download_eggnog_data.py and follow '
                                                                       'prompts. The Diamond AND HMM databases are '
                                                                       'required.')
        exit(OSError)
