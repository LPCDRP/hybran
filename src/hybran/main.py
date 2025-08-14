import os
import argparse
import logging
import time
import glob
import tempfile
import atexit
import shutil
import warnings

from . import (
    verifyInstallations,
    fileManager,
    extractor,
    onegene,
    run,
    annomerge,
    converter,
    config,
    designator,
    standardize,
    compare,
    defuse,
    __version__,
)
from .argparse import DefaultSubcommandArgumentParser


warnings.filterwarnings(
    "ignore",
    message=r"Partial codon",
)

def percentage(string):
    num = float(string)
    if not (num <= 100 and num >= 0):
        raise argparse.ArgumentTypeError("Not a percentage.")
    return num

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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    head_parser.set_default_subparser('annotate')


    # alternative commands
    stdize = subparsers.add_parser(
        'standardize',
        help='Apply standard naming conventions to Hybran output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stdize.set_defaults(func=standardize.main)
    onegenecmd = subparsers.add_parser(
        'onegene',
        help='Unify names of gene duplicates.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    onegenecmd.set_defaults(func=onegene.main)
    comparecmd = subparsers.add_parser(
        'compare',
        help='Compare two annotations of the same genome.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    comparecmd.set_defaults(func=compare.main)
    defusecmd = subparsers.add_parser(
        'defuse',
        help='Separate gene fusion annotations into single gene annotations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    defusecmd.set_defaults(func=defuse.main)


    #
    # hybran standardize
    #
    stdize.add_argument(
        'annotations',
        help=(
            "Directory, space-separated list of GBKs, or a FOFN containing all annotated genomes. "
            "If you pass a hybran output directory here, no other arguments will be required. "
            "Otherwise, you will need to provide the location of the unifications file via -u/--unifications-file."
        ),
        nargs='+'
    )
    stdize.add_argument(
        '-p', '--orf-prefix',
        type=str,
        help='prefix for generic gene names (*not* locus tags)',
        default='HYBRA',
    )
    stdize.add_argument(
        '-o', '--output',
        help='Directory to output all new annotation files.',
        default='.',
    )
    stdize.add_argument(
        '-u', '--unifications-file',
        help="reference annotation's unifications.tsv file produced by hybran onegene.",
    )
    stdize.add_argument(
        '-r', '--ref-names-only',
        action='store_true',
        help="Do not use gene names supplied by the ab initio caller.",
    )

    #
    # hybran onegene
    #
    onegenecmd.add_argument(
        'annotations',
        help="Directory, space-separated list of GBKs, or a FOFN containing all annotated genomes.",
        nargs='+'
    )
    onegenecmd.add_argument(
        '-p', '--orf-prefix',
        type=str,
        help=(
            "prefix for unifying gene names (*not* locus tags). "
            "Such names will be applied to all sets of highly conserved genes if they don't already have a name or if they have discrepant names. "
            "Whatever you pass here will be sandwiched by REF and X. (i.e., the default HYBRA will be transformed into REFHYBRAX and then used)."
        ),
        default='HYBRA',
    )
    onegenecmd.add_argument(
        '-o', '--output',
        help='Directory to output all new annotation files.',
        default='.',
    )
    onegenecmd.add_argument(
        '-i', '--identity-threshold',
        required=False,
        type=int,
        help='Percent sequence identity threshold to use during CD-HIT clustering and BLASTP to call duplications',
        default=99,
    )
    onegenecmd.add_argument(
        '-c', '--coverage-threshold',
        required=False,
        type=int,
        help='Percent alignment coverage threshold to use during CD-HIT clustering and BLASTP to call duplications',
        default=99,
    )
    onegenecmd.add_argument(
        '-t', '--first-reference',
        required=False,
        dest='first_gbk',
        help="Reference name or file name whose locus tags should be used as unified names for conserved copies in the others."
        " Default is the annotation with the most named CDSs. If you specify a file here that is not in your input list, it will be added."
    )


    #
    # hybran compare
    #
    comparecmd.add_argument(
        'annotations',
        help='The two annotation files to compare, in genbank format.',
        nargs=2
    )
    comparecmd.add_argument(
        '-o', '--outdir',
        help='Directory to output the results of the comparison.',
        default='.',
    )

    #
    # hybran defuse
    #
    defusecmd.add_argument(
        'annotations_dir',
        help=(
            "Results directory from the Hybran run whose gene fusions you wish to defuse."
        ),
    )
    defusecmd.add_argument(
        '-o', '--output',
        help='Directory to output all new annotation files.',
        default='.',
    )

    #
    # hybran annotate
    #
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    bbh_parser = parser.add_argument_group(
        'Bidirectional BLAST Hit (BBH) Settings',
        description=(
            'BBHs are used in matching ab initio predicted genes to reference genes,'
            ' as well as to resolve members of MCL clusters containing multiple different reference genes.'
        ),
    )
    ratt_params = parser.add_argument_group('RATT Options.\n(See http://ratt.sourceforge.net/documentation.html and\n https://github.com/ThomasDOtto/ratt for more details)')
    prokka_params = parser.add_argument_group('Prokka Options.\n(See https://github.com/tseemann/prokka for more details)')
    required.add_argument('-g', '--genomes', help='Directory, a space-separated list of FASTAs, or a FOFN '
                                                      'containing all genomes desired to be annotated. '
                                                      'FASTA format required.',
                          required=True,
                          nargs='+')
    required.add_argument('-r', '--references', help='Directory, a space-separated list of GBKs, or a FOFN '
                                                     'containing Genbank files of reference '
                                                     'annotations to transfer.',
                          required=True,
                          nargs='+')
    # This flag is set to match the way pgap takes it.
    optional.add_argument('-s', '--organism',
                          help="genus only or binomial name",
                          required=False,
                          )
    optional.add_argument('-e', '--eggnog-databases', help='Directory of the eggnog databases downloaded using '
                                                           'download_eggnog_data.py -y bactNOG. Full path only',
                          dest='database_dir',
                          required=False)
    optional.add_argument('-t', '--first-reference', required=False, dest='first_gbk',
                          help="Reference name or file name whose locus tags should be used as unified names for conserved copies in the others."
                          " Default is the annotation with the most named CDSs. If you specify a file here that is not in your input list, it will be added.")
    optional.add_argument(
        '-I', '--mcl-inflation',
        help=(
            "MCL inflation value. "
            "Higher value results in more fine-grained clusters (fewer genes in common). "
            "See <https://micans.org/mcl/man/mcl.html#opt-I> for details."
        ),
        type=float,
        default=config.cnf.mcl_inflation,
    )
    optional.add_argument('-i', '--identity-threshold', required=False, type=int,
                          help='Percent sequence identity threshold to use during CD-HIT clustering and BLASTP',
                          default=99)
    optional.add_argument('-c', '--coverage-threshold', required=False, type=int,
                          help='Percent sequence coverage threshold to use during CD-HIT clustering and BLASTP',
                          default=99)
    optional.add_argument('-o', '--output', help='Directory to output all new annotation files.',
                          default='.')
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
                          default='HYBRA')
    optional.add_argument('--dedupe-references',
                          action='store_true',
                          help=(
                              'Identify duplicate genes in the reference annotations and assign one name to all copies.'
                              'This option is deprecated and has no effect.'
                              'If name unification is not desired, consider running `hybran standardize` afterwards.'
                          ))

    bbh_parser.add_argument(
        '--bbh-min-coverage',
        help="Minimum percent query and subject alignment coverage for candidate BBHs.",
        default=config.cnf.blast.min_coverage,
        type=percentage,
    )
    bbh_parser.add_argument(
        '--bbh-min-bitscore',
        help="Minimum BLAST bitscore for candidate BBHs.",
        default=config.cnf.blast.min_bitscore,
        type=float,
    )

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
                                 'Strain.global',
                                 'Strain.global.Repetitive',
                                 'Species',
                                 'Species.Repetitive',
                                 'Species.global',
                                 'Species.global.Repetitive',
                                 'Pacbio',
                                 'Pacbio.Repetitive',
                                 'PacbioG',
                                 'Falciparum',
                                 'Falciparum.Repetitive',
                                 'Multiple',
                                 'Free',
                             ],
                             default='Strain',
                             help=(
                                 'Presets for nucmer alignment settings to determine synteny.'
                                 "Automatically set to 'Multiple' when multiple references are provided unless 'Free' is specified."
                             )
    )
    ratt_params.add_argument('--ratt-splice-sites',
                             type=str,
                             default='XX..XX',
                             help='splice donor and acceptor sequences. example: GT..AG')
    ratt_params.add_argument('--ratt-correct-splice',
                             action=argparse.BooleanOptionalAction,
                             default=False,
                             help='whether RATT should attempt splice site corrections')
    ratt_params.add_argument('--ratt-correct-pseudogenes',
                             action=argparse.BooleanOptionalAction,
                             default=False,
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
    prokka_params.add_argument('--genus', help='Genus name. Deprecated -- use hybran -s/--organism instead')
    prokka_params.add_argument('--species', help='Species name. Deprecated -- use hybran -s/--organism instead')
    prokka_params.add_argument('--strain', help='Strain name. Deprecated -- use hybran -s/--organism instead')
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
    config.cnf.mcl_inflation = args.mcl_inflation
    config.cnf.blast.min_coverage = args.bbh_min_coverage
    config.cnf.blast.min_bitscore = args.bbh_min_bitscore

    designator.generic_orf_prefix[0]=args.orf_prefix
    designator.ref_orf_prefix[0] = f"REF{args.orf_prefix}X"

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

    # Convert all input paths to full path if not given as full path
    args.genomes = fileManager.file_list(args.genomes, file_type="fasta")
    args.references = fileManager.file_list(args.references, file_type="genbank")
    args.output = os.path.abspath(args.output)

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit(f"Could not create directory {args.output}")

    # Setting up logging
    start_time = time.time()
    if not args.quiet:
        if args.verbose:
            config.cnf.verbosity = 2
            loglevel = logging.DEBUG
        else:
            config.cnf.verbosity = 1
            loglevel = logging.INFO
    else:
        config.cnf.verbosity = 0
        loglevel = logging.WARNING
    logfile = os.path.join(args.output, "hybran.log")
    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(),
        ],
        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
        datefmt='%H:%M:%S %m-%d',
    )
    rootlogger = logging.getLogger()
    [filehandler, streamhandler] = rootlogger.handlers
    streamhandler.setLevel(loglevel)

    logger = logging.getLogger('Hybran')

    if not args.organism:
        organism = "Genus species"
        strain = ""
    else:
        organism = args.organism
        strain = ""
        name_components = organism.split()
        if len(name_components) > 2:
            organism = name_components[0:2]
            strain = ' '.join(name_components[2:])

    if args.dedupe_references:
        logger.warning("""The --dedupe-references option is deprecated, as name unification happens constitutively.
This option is scheduled for removal, so please update your invocation for the future.
"""
        )

    # Moving into the desired annotation directory
    os.chdir(args.output)

    # Setting up references for RATT, as well as versions in GFF format used later
    if not os.path.isdir('unified-refs'):
        try:
            os.mkdir('unified-refs')
        except:
            sys.exit("Could not create directory: unified-refs ")
    if args.first_gbk is not None and os.path.isfile(args.first_gbk):
        if os.path.abspath(args.first_gbk) not in [
                os.path.abspath(_) for _ in args.references
        ]:
            args.references.append(args.first_gbk)
            args.first_gbk = os.path.basename(os.path.splitext(args.first_gbk)[0])
    deduped_refs = [os.path.abspath(os.path.join('unified-refs',os.path.basename(_))) for _ in args.references]
    ref_cds = os.path.abspath(os.path.join('unified-refs', 'unique_ref_cdss.faa'))
    if not all([os.path.isfile(_) for _ in deduped_refs]):
        ref_cds, args.references = onegene.unify(
            args.references,
            outdir='unified-refs',
            tmpdir=hybran_tmp_dir,
            main_ref=args.first_gbk,
        )
    else:
        args.references = deduped_refs
    refdir, embl_dir, embls = fileManager.prepare_references(args.references)
    ref_gbks = fileManager.file_list([refdir], "genbank")
    if len(ref_gbks) > 1 and args.ratt_transfer_type not in ["Multiple", "Free"]:
        logger.info(f"Multiple reference annotations found. Setting RATT transfer type to 'Multiple' accordingly...")
        args.ratt_transfer_type = "Multiple"
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
    genetic_code = extractor.get_genetic_code(ref_gbks[0])
    config.set_genetic_code(genetic_code)
    logger.info('Using genetic code ' + str(genetic_code) + ' as detected in reference annotation.')


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
                                organism=organism,
                                strain=strain,
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
                                  organism=organism,
                                  strain=strain,
                                  genome=genome,
                                  annotation_fp=os.getcwd() + '/',
                                  ref_proteins_fasta=ref_cds,
                                  ref_gbk_list=ref_gbks,
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
        if args.database_dir and 'eggnog_seqs.fasta' in os.listdir(hybran_tmp_dir):
            ref_tax_ids = [extractor.get_taxonomy_id(_) for _ in ref_gbks]
            if any(ref_tax_ids):
                ref_features_by_record = {}
                for ref_gbk_file in ref_gbks:
                    refname = os.path.splitext(os.path.basename(ref_gbk_file))[0]
                    for record in SeqIO.parse(ref_gbk_file, "genbank"):
                        ref_features_by_record[f"{refname}.{record.id}"] = record.features
                # TODO: extractor.gene_dict should have another level by refname
                ref_gene_dict = extractor.gene_dict(ref_features_by_record)
                [ref_gene_dict.update(extractor.gene_dict(_)) for _ in ref_gbks]
                run.eggnog_mapper(
                    script_dir=script_dir,
                    nproc=args.nproc,
                    emapper_loc=args.database_dir,
                    ref_tax_ids = [_ for _ in ref_tax_ids if _ is not None],
                    ref_gene_dict = ref_gene_dict,
                    temp_dir=hybran_tmp_dir
                )
            else:
                logger.info('No reference taxonomy IDs detected. Skipping eggNOG mapping...')
        else:
            logger.info('No genes to be annotated with eggnog, continuing')

        logger.info('Finalizing annotations...')
        for gbk in glob.glob(os.path.join(args.output,'*.gbk')):
            designator.create_gene_entries(gbk)
            converter.convert_gbk_to_gff(gbk)

    logger.info('Finished. Annotated ' + str(genome_count) + ' genomes. Genbank and GFF are located in ' + args.output)
    logger.info('Time elapsed: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes\n')

if __name__ == '__main__':
    cmds()
