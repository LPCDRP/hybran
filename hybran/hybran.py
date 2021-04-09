import os
import argparse
import logging
import time
import glob
import tempfile
import atexit
import shutil

from . import \
    verifyInstallations, \
    fileManager, \
    firstReference, \
    run, \
    annomerge, \
    converter, \
    config


__version__ = '1.1.1'


def cmds():
    """
    argparse parse input provided by the user

    :return: argparse.parse_args() object
    """
    parser = argparse.ArgumentParser(description='Hybran: a pipeline to annotate Mycobacterium '
                                                 'tuberculosis de novo assembled genomes. Annotation of other species '
                                                 'or mixing species within an annotation run is NOT recommended.'
                                                 '\n\nPlease cite: [manuscript submitted]')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    required.add_argument('-g', '--genomes', help='Directory, a space-separated list of FASTAs, or a FOFN '
                                                      'containing all genomes desired to be annotated. '
                                                      'FASTA format required.',
                          required=True,
                          nargs='+')
    required.add_argument('-r', '--references', help='Directory containing Genbank files of reference '
                                                     'annotations to transfer. Only the first annotation will '
                                                     'be used as the reference database in the Prokka reference '
                                                     'step.',
                          required=True)
    required.add_argument('-e', '--eggnog-databases', help='Directory of the eggnog databases downloaded using '
                                                           'download_eggnog_data.py -y bactNOG. Full path only',
                          dest='database_dir',
                          required=True)
    optional.add_argument('-t', '--first-reference', required=False, dest='first_gbk',
                          help='Reference to use as the reference database for Prokka. Must exist in --references dir.'
                               ' Default is the first reference annotation (Genbank) in -r/--references.')
    optional.add_argument('-i', '--identity-threshold', required=False, type=int,
                          help='Percent sequence identity threshold to use during CD-HIT clustering and BLASTP. Default is 95',
                          default=95)
    optional.add_argument('-c', '--coverage-threshold', required=False, type=int,
                          help='Percent sequence coverage threshold to use during CD-HIT clustering and BLASTP. Default is 95',
                          default=95)
    optional.add_argument('-o', '--output', help='Directory to output all new annotation files. Default is the '
                                                 '-r/--references directory. Full path only')
    optional.add_argument('-n', '--nproc', help='Number of processors/CPUs to use. Default is 1',
                          default='1')
    optional.add_argument('-d', '--debug', action='store_true',
                          help="Don't delete temporary files created by Hybran.")
    optional.add_argument('-f', '--force', action='store_true',
                          help='Force overwrite intermediate files (does not overwrite annotation files already '
                               'annotated using hybran.')
    logging_level = optional.add_mutually_exclusive_group()
    logging_level.add_argument('--verbose', action='store_true', help='Verbose output')
    logging_level.add_argument('-q', '--quiet', action='store_true', help='No logging output when flagged')
    optional.add_argument('-v', '--version', help='Print version and exit',
                          action='version',
                          version=__version__)
    return parser.parse_args()


def main():
    """
    Hybran: a pipeline to annotate Mycobacterium
    tuberculosis de novo assembled genomes. Annotation of other species
    or mixing species within an annotation run is NOT recommended.

    :return: None
    """
    global args
    args = cmds()
    # Obtaining the absolute path to all scripts
    script_dir = os.path.abspath(os.path.dirname(__file__))

    # Setting up the Hybran temporary directory
    config.init()
    hybran_tmp_dir = config.hybran_tmp_dir

    # Cleanup the temporary files directory and its contents at exit unless
    # --debug is set
    if not args.debug:
        atexit.register(shutil.rmtree, path=hybran_tmp_dir)

    # Check that the identity threshold is valid
    if not (args.identity_threshold <= 100 and args.identity_threshold >= 0):
        print("error: invalid value for --identity-threshold. Must be between 0 and 100.")
        exit(10)

    # Check that the coverage threshold is valid
    if not (args.coverage_threshold <= 100 and args.coverage_threshold >= 0):
        print("error: invalid value for --coverage-threshold. Must be between 0 and 100.")
        exit(10)

    # Confirming all installations are valid
    verifyInstallations.verify_installations(args.database_dir)

    # Setting up logging
    start_time = time.time()
    print('\n\t\t\tPlease cite:\n\t\t\t[manuscript submitted]\n\n')
    if not args.quiet:
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
                                datefmt='%H:%M:%S %m-%d')
        else:
            logging.basicConfig(level=logging.INFO,
                                format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
                                datefmt='%H:%M:%S %m-%d')
    logger = logging.getLogger('Hybran')
    if not args.output:
        args.output = args.references
    cwd = os.getcwd() + '/'

    # Convert all input paths to full path if not given as full path
    if os.path.isdir(args.genomes[0]):
        genomes = os.listdir(fileManager.full_path(args.genomes[0]))
        args.genomes = [args.genomes[0] + g for g in genomes]
    elif len(args.genomes) > 1:
        args.genomes = [fileManager.full_path(g) for g in args.genomes]
    elif args.genomes[0].endswith('.fasta') or args.genomes[0].endswith('.fna') or args.genomes[0].endswith('.fa'):
        args.genomes = [fileManager.full_path(args.genomes[0])]
    else:
        genomes = []
        with open(args.genomes[0], 'r') as genome:
            for line in genome:
                genomes.append(line.rstrip('\n'))
        args.genomes = genomes
    args.references = fileManager.full_path(args.references)
    args.output = fileManager.full_path(args.output)

    # Moving into the desired annotation directory
    os.chdir(args.output)

    # Setting up references for RATT, as well as versions in GFF format used later
    refdir, embl_dir, embls = fileManager.prepare_references(args)
    intermediate_dirs = ['clustering/', 'eggnog-mapper-annotations/', 'prodigal-test/', refdir] + \
                        [d for d in glob.glob('emappertmp*/')]
    intermediate_files = ['reference_prodigal_proteome.faa', 'ref.fasta', 'eggnog_seqs.fasta', 'ref_proteome.fasta']
    if all(intermediate_files + intermediate_dirs) in os.listdir(hybran_tmp_dir) and \
            not args.force:
        print('It seems hybran has been executed before. Please use the '
              '-f/--force flag if you would like to ' \
              'overwrite intermediate files ' \
              '(does not overwrite Genbank/GFF files created in a previous run of hybran).')
        exit()
    # Getting first reference information
    if not args.first_gbk:
        first_reference_embl = embls[0]
        first_reference_gbk = os.path.splitext(first_reference_embl)[0] + '.gbk'
    else:
        first_reference_gbk = os.path.join(refdir, args.first_gbk)
        first_reference_embl = os.path.splitext(args.first_gbk)[0] + '.embl'
    ref_cds, ref_genome = firstReference.get_first_reference_proteome(first_reference_gbk)

    # Calling all steps for Hybran
    genome_count = 0
    genomes_annotate = []
    all_genomes = []
    genomes = []
    for genome in args.genomes:
        if genome.endswith('.fasta') and not genome.startswith('ref'):
            filename = os.path.splitext(genome)[0]
            samplename = os.path.basename(filename)
            if samplename + '.gbk' not in args.output:
                genome_count += 1
                genomes.append(filename)
                run.ratt_prokka(ref_dir=embl_dir,
                                fasta=genome,
                                ref_cds=ref_cds,
                                script_dir=script_dir,
                                cpus=args.nproc)
                if samplename + '.gbk' not in os.listdir(os.getcwd()):
                    logger.info('Merging RATT and Prokka annotations for ' + samplename)
                    annomerge.run(isolate_id=samplename,
                                  annotation_fp=os.getcwd() + '/',
                                  ref_proteins_fasta=ref_cds,
                                  ref_embl_fp=first_reference_embl,
                                  reference_genome=ref_genome,
                                  script_directory=script_dir,
                                  seq_ident=args.identity_threshold,
                                  seq_covg=args.coverage_threshold)
                genomes_annotate.append(os.path.join(args.output, samplename + '.gff'))
            else:
                logger.info(genome + ' as already been annotated.')
            if samplename + '.gff' not in os.listdir(args.output):
                converter.convert_gbk_to_gff(os.path.join(args.output, samplename + '.gbk'))
    if all([os.path.isfile(os.path.join(args.output, os.path.basename(g) + '.gff')) for g in genomes]):
        all_genomes += [refdir + i for i in os.listdir(refdir) if i.endswith('.gff')] + genomes_annotate
        run.clustering(all_genomes=list(set(sorted(all_genomes))),
                       target_genomes=genomes_annotate,
                       nproc=args.nproc,
                       seq_ident=args.identity_threshold,
                       seq_covg=args.coverage_threshold)
        if 'eggnog_seqs.fasta' in os.listdir(hybran_tmp_dir):
            run.eggnog_mapper(script_dir=script_dir,
                              nproc=args.nproc,
                              emapper_loc=args.database_dir,
                              temp_dir=hybran_tmp_dir)
        else:
            logger.info('No genes to be annotated with eggnog, continuing')

        for gbk in os.listdir(args.output):
            if gbk.endswith('.gbk'):
                converter.convert_gbk_to_gff(gbk)

    logger.info('Finished. Annotated ' + str(genome_count) + ' genomes. Genbank and GFF are located in ' + args.output)
    logger.info('Time elapsed: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes\n')
    print('Thank you for using Hybran. We hope to see you again!')


if __name__ == '__main__':
    main()
