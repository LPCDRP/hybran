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
    extractor, \
    refManager, \
    run, \
    annomerge, \
    converter, \
    config, \
    designator, \
    standardize, \
    __version__

from .argparse import DefaultSubcommandArgumentParser


def cmds():
    """
    argparse parse input provided by the user and call appropriate command
    """
    citation = (
        "Elghraoui, A.; Gunasekaran, D.; Ramirez-Busby, S. M.; Bishop, E.; Valafar, F. "
        "Hybran: Hybrid Reference Transfer and Ab Initio Prokaryotic Genome Annotation. "
        "bioRxiv November 10, 2022, p 2022.11.09.515824. "
        "<https://doi.org/10.1101/2022.11.09.515824>"
    )

    head_parser = DefaultSubcommandArgumentParser(
        description='Hybran: hybrid reference-based and ab initio prokaryotic genomic annotation.',
        epilog=citation,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = head_parser.add_subparsers(
        dest='subparser_name',
    )
    # We want the main annotation command to run as default, but also provide alternative commands.
    #
    # Standard argparse seems to handle this, with this `required` option to `add_subparsers`,
    # but it's not actually what we want.
    #  > required - Whether or not a subcommand must be provided, by default False (added in 3.7)
    # subparsers = parser.add_subparsers(
    #     title="Alternative subcommands",
    #     required=False,
    # )

    # main command "annotate"
    parser = subparsers.add_parser(
        'annotate',
        help='Run the Hybran annotation pipeline. (default)',
    )
    head_parser.set_default_subparser('annotate')


    # alternative commands
    stdize = subparsers.add_parser(
        'standardize',
        help='Apply standard naming conventions to Hybran output.'
    )
    stdize.set_defaults(func=standardize.main)


    #
    # hybran standardize
    #
    stdize.add_argument(
        'annotations',
        help="Directory, space-separated list of GBKs, or a FOFN containing all annotated genomes.",
        nargs='+'
    )
    stdize.add_argument(
        '-p', '--orf-prefix',
        type=str,
        help='prefix for generic gene names (*not* locus tags)',
        default='ORF',
    )
    stdize.add_argument(
        '-o', '--output',
        help='Directory to output all new annotation files.',
        default='.',
    )

    #
    # hybran annotate
    #
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    ratt_params = parser.add_argument_group('RATT Options.\n(See http://ratt.sourceforge.net/documentation.html and\n https://github.com/ThomasDOtto/ratt/blob/master/ratt.1.md for more details)')
    prokka_params = parser.add_argument_group('Prokka Options.\n(See https://github.com/tseemann/prokka for more details)')
    required.add_argument('-g', '--genomes', help='Directory, a space-separated list of FASTAs, or a FOFN '
                                                      'containing all genomes desired to be annotated. '
                                                      'FASTA format required.',
                          required=True,
                          nargs='+')
    required.add_argument('-r', '--references', help='Directory, a space-separated list of GBKs, or a FOFN '
                                                     'containing Genbank files of reference '
                                                     'annotations to transfer. Only the first annotation will '
                                                     'be used as the reference database in the Prokka reference '
                                                     'step.',
                          required=True,
                          nargs='+')
    optional.add_argument('-e', '--eggnog-databases', help='Directory of the eggnog databases downloaded using '
                                                           'download_eggnog_data.py -y bactNOG. Full path only',
                          dest='database_dir',
                          required=False)
    optional.add_argument('--dedupe-references',
                          action='store_true',
                          help='Identify duplicate genes in the reference annotations and assign one name to all copies.')
    optional.add_argument('-t', '--first-reference', required=False, dest='first_gbk',
                          help='Reference to use as the reference database for Prokka. Must exist in --references dir.'
                               ' Default is the first reference annotation (Genbank) in -r/--references.')
    optional.add_argument('-i', '--identity-threshold', required=False, type=int,
                          help='Percent sequence identity threshold to use during CD-HIT clustering and BLASTP',
                          default=95)
    optional.add_argument('-c', '--coverage-threshold', required=False, type=int,
                          help='Percent sequence coverage threshold to use during CD-HIT clustering and BLASTP',
                          default=95)
    optional.add_argument('-o', '--output', help='Directory to output all new annotation files. Default is the '
                                                 '-r/--references directory. Full path only')
    optional.add_argument('-n', '--nproc', help='Number of processors/CPUs to use',
                          type=int,
                          default=1)
    optional.add_argument('-d', '--debug', action='store_true',
                          help="Don't delete temporary files created by Hybran.")
    optional.add_argument('-f', '--force', action='store_true',
                          help='Force overwrite intermediate files (does not overwrite annotation files already '
                               'annotated using hybran.')
    optional.add_argument('--filter-ratt', action='store_true',
                          help='Enforce identity/coverage thresholds on RATT-transferred annotations.')
    optional.add_argument('-p', '--orf-prefix',
                          type=str,
                          help='prefix for generic gene names (*not* locus tags)',
                          default='ORF')
    logging_level = optional.add_mutually_exclusive_group()
    logging_level.add_argument('--verbose', action='store_true', help='Verbose output')
    logging_level.add_argument('-q', '--quiet', action='store_true', help='No logging output when flagged')
    optional.add_argument('-v', '--version', help='Print version and exit',
                          action='version',
                          version=__version__)

    ratt_params.add_argument('--ratt-transfer-type',
                             choices=[
                                 'Assembly',
                                 'Assembly.Repetitive',
                                 'Strain',
                                 'Strain.Repetitive',
                                 'Strain.Global',
                                 'Species',
                                 'Species.Repetitive',
                                 'Species.Global',
                                 'Multiple',
                                 'Free',
                             ],
                             default='Strain',
                             help='presets for nucmer alignment settings to determine synteny.')
    ratt_params.add_argument('--ratt-splice-sites',
                             type=str,
                             default='XX..XX',
                             help='splice donor and acceptor sequences. example: GT..AG')
    ratt_params.add_argument('--ratt-correct-splice',
                             action = 'store_true',
# TODO: need python 3.9 for this
#                             action = argparse.BooleanOptionalAction,
                             help='whether RATT should attempt splice site corrections')
    ratt_params.add_argument('--ratt-correct-pseudogenes',
                             action = 'store_true',
# TODO: need python 3.9 for this
#                             action = argparse.BooleanOptionalAction,
                             help='whether RATT should attempt correction of reference pseudogenes in your samples')


    prokka_params.add_argument('--kingdom',
                               choices=[
                                   'Archaea',
                                   'Bacteria',
                                   'Mitochondria',
                                   'Viruses',
                               ],
                               default='Bacteria',
                               help='Determines which UniProtKB databases Prokka searches against.')
    prokka_params.add_argument('--genus', help='Genus name')
    prokka_params.add_argument('--species', help='Species name')
    prokka_params.add_argument('--strain', help='Strain name')
    prokka_params.add_argument('--plasmid', help='Plasmid name or identifier')
    prokka_params.add_argument('--gram',
                               choices=[
                                   '+',
                                   'pos',
                                   '-',
                                   'neg',
                               ],
                               help='Gram')
    prokka_params.add_argument('--prodigaltf',
                               help="Prodigal training file")
    prokka_params.add_argument('--hmms',
                               help='Trusted HMM to first annotate from')
    prokka_params.add_argument('--metagenome',
                               action='store_true',
                               help="Improve gene predictions for highly fragmented genomes")
    prokka_params.add_argument('--evalue',
                               type=float,
                               default=1e-9,
                               help='Similarity e-value cut-off')

    arguments = head_parser.parse_args()

    if arguments.subparser_name != 'annotate':
        arguments.func(arguments)
        return

    # Turn Prokka arguments back into an option string that we can pass directly to Prokka,
    # rather than having to repopulate --kingdom args.kingdom --gram args.gram ...
    #
    # special thanks: https://stackoverflow.com/a/31520622
    prokka_passthrough = ' '.join(
        [' '.join([ g.option_strings[0], str(vars(arguments)[g.dest]) ])
         # action_groups[-1] is the last defined argument_group, i.e. prokka
         for g in parser._action_groups[-1]._group_actions if vars(arguments)[g.dest]]
    ).replace("True","")

    main(arguments, prokka_passthrough)


def main(args, prokka_args):
    """
    Hybran: a pipeline to annotate Mycobacterium
    tuberculosis de novo assembled genomes. Annotation of other species
    or mixing species within an annotation run is NOT recommended.

    :return: None
    """
    # Obtaining the absolute path to all scripts
    script_dir = os.path.abspath(os.path.dirname(__file__))

    # Setting up the Hybran temporary directory
    config.init()
    hybran_tmp_dir = config.hybran_tmp_dir

    designator.generic_orf_prefix[0]=args.orf_prefix

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
    if args.database_dir:
        verifyInstallations.verify_installations(args.database_dir)

    # Setting up logging
    start_time = time.time()
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
    args.genomes = fileManager.file_list(args.genomes, file_type="fasta")
    args.references = fileManager.file_list(args.references, file_type="genbank")
    args.output = os.path.abspath(args.output)

    # Moving into the desired annotation directory
    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit("Could not create directory " + args.output)
    os.chdir(args.output)

    # Setting up references for RATT, as well as versions in GFF format used later
    if args.dedupe_references:
        if not os.path.isdir('deduped-refs'):
            try:
                os.mkdir('deduped-refs')
            except:
                sys.exit("Could not create directory: deduped-refs ")
        args.references = refManager.dedupe(args.references, outdir='deduped-refs', tmpdir=hybran_tmp_dir)
    refdir, embl_dir, embls = fileManager.prepare_references(args.references)
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
        first_reference_gbk = glob.glob(os.path.join(refdir, '*.gbk'))[0]
    else:
        first_reference_gbk = os.path.join(refdir, args.first_gbk)
    ref_cds     = os.path.join(hybran_tmp_dir, 'ref_proteome.fasta')
    ref_genome  = os.path.join(hybran_tmp_dir, 'ref.fasta')
    genetic_code = extractor.get_genetic_code(first_reference_gbk)
    config.set_genetic_code(genetic_code)
    logger.info('Using genetic code ' + str(genetic_code) + ' as detected in reference annotation.')
    logger.info('Creating a reference proteome FASTA for Prokka from ' + first_reference_gbk)
    extractor.fastaFromGbk(
        genbank = first_reference_gbk,
        out_cds = ref_cds,
        out_genome = ref_genome,
        describe = extractor.prokka_faa,
    )


    # Configure RATT start/stop codons and other settings,
    # but only if the user hasn't defined their own
    if 'RATT_CONFIG' not in os.environ:
        os.environ['RATT_CONFIG'] = os.path.join(hybran_tmp_dir,'RATT.config')
        config.ratt(
            os.environ['RATT_CONFIG'],
            genetic_code,
            splice = args.ratt_splice_sites,
            correct_splice = args.ratt_correct_splice,
            correct_pseudo = args.ratt_correct_pseudogenes,
        )


    # Calling all steps for Hybran
    genome_count = 0
    genomes_annotate = []
    all_genomes = []
    genomes = []
    for genome in args.genomes:
        if any([genome.endswith('.' + ext) for ext in fileManager.exts['fasta']]) and genome != 'ref.fasta':
            filename = os.path.splitext(genome)[0]
            samplename = os.path.basename(filename)
            annomerge_gbk = os.path.join(args.output, samplename, 'annomerge', samplename + '.gbk')
            annomerge_gff = os.path.join(args.output, samplename, 'annomerge', samplename + '.gff')
            if samplename + '.gbk' not in args.output:
                genome_count += 1
                genomes.append(filename)
                run.ratt_prokka(ref_dir=embl_dir,
                                fasta=genome,
                                ref_cds=ref_cds,
                                gcode=genetic_code,
                                ratt_ttype = args.ratt_transfer_type,
                                prokka_extra_args = prokka_args,
                                script_dir=script_dir,
                                cpus=args.nproc,
                                qcov=args.coverage_threshold)
                if not os.path.isfile(annomerge_gbk):
                    logger.info('Merging RATT and Prokka annotations for ' + samplename)
                    annomerge.run(isolate_id=samplename,
                                  contigs=extractor.get_contig_names(genome),
                                  annotation_fp=os.getcwd() + '/',
                                  ref_proteins_fasta=ref_cds,
                                  ref_gbk_fp=first_reference_gbk,
                                  reference_genome=ref_genome,
                                  script_directory=script_dir,
                                  seq_ident=args.identity_threshold,
                                  seq_covg=args.coverage_threshold,
                                  ratt_enforce_thresholds=args.filter_ratt,
                                  nproc=args.nproc,
                    )
                    converter.convert_gbk_to_gff(annomerge_gbk)
                genomes_annotate.append(os.path.join(args.output, samplename + '.gff'))
            else:
                logger.info(genome + ' as already been annotated.')
            if samplename + '.gff' not in os.listdir(args.output):
                shutil.copy(annomerge_gbk, args.output)
                shutil.copy(annomerge_gff, args.output)
    if all([os.path.isfile(os.path.join(args.output, os.path.basename(g) + '.gff')) for g in genomes]):
        all_genomes += [refdir + i for i in os.listdir(refdir) if i.endswith('.gff')] + genomes_annotate
        run.clustering(all_genomes=list(set(sorted(all_genomes))),
                       target_genomes=genomes_annotate,
                       nproc=args.nproc,
                       seq_ident=args.identity_threshold,
                       seq_covg=args.coverage_threshold)
        tax_id = extractor.get_taxonomy_id(first_reference_gbk)
        if args.database_dir and tax_id and 'eggnog_seqs.fasta' in os.listdir(hybran_tmp_dir):
            ref_gene_dict = extractor.gene_dict(first_reference_gbk)
            run.eggnog_mapper(script_dir=script_dir,
                              nproc=args.nproc,
                              emapper_loc=args.database_dir,
                              ref_tax_id = tax_id,
                              ref_gene_dict = ref_gene_dict,
                              temp_dir=hybran_tmp_dir)
        else:
            logger.info('No genes to be annotated with eggnog, continuing')

        logger.info('Assigning locus tags')
        for gbk in glob.glob(os.path.join(args.output,'*.gbk')):
            isolate_id = os.path.basename(os.path.splitext(gbk)[0])
            designator.assign_locus_tags(gbk, prefix=isolate_id)
            designator.create_gene_entries(gbk)
            converter.convert_gbk_to_gff(gbk)

    logger.info('Finished. Annotated ' + str(genome_count) + ' genomes. Genbank and GFF are located in ' + args.output)
    logger.info('Time elapsed: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes\n')
    print('Thank you for using Hybran. We hope to see you again!')


if __name__ == '__main__':
    cmds()
