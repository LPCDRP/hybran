import atexit
from collections import defaultdict
import os
import shutil
import sys

from Bio import SeqIO

from . import config, fileManager
from . import CDHIT, designator, extractor, parseClustering


def main(args):
    config.init()
    atexit.register(shutil.rmtree, path=config.hybran_tmp_dir)

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit("Could not create directory " + args.output)

    # Check that the identity threshold is valid
    if not (args.identity_threshold <= 100 and args.identity_threshold >= 0):
        sys.exit("error: invalid value for --identity-threshold. Must be between 0 and 100.")

    # Check that the coverage threshold is valid
    if not (args.coverage_threshold <= 100 and args.coverage_threshold >= 0):
        sys.exit("error: invalid value for --coverage-threshold. Must be between 0 and 100.")


    annotations = fileManager.file_list(args.annotations, file_type="genbank")
    main_ref = args.first_gbk
    if main_ref is None or main_ref in [os.path.basename(os.path.splitext(_)[0]) for _ in annotations]:
        pass
    elif os.path.isfile(os.path.abspath(main_ref)):
        main_ref_file = os.path.abspath(main_ref)
        if main_ref_file not in annotations:
            annotations.append(main_ref_file)
            main_ref = os.path.splitext(os.path.basename(main_ref))[0]
    else:
        sys.exit(f"error: could not find primary reference annotation {main_ref}")

    unify(
        annotations=annotations,
        outdir=args.output,
        tmpdir=config.hybran_tmp_dir,
        seq_ident=args.identity_threshold,
        seq_covg=args.coverage_threshold,
        main_ref=main_ref,
    )

def unify(annotations, outdir, tmpdir, seq_ident=99, seq_covg=99, main_ref=None):
    """"
    Identify conserved/duplicated genes among the reference annotation and assign all instances a
    single name. If instances had different actual gene names, these original names are
    stored in the gene_synonym field.

    :param annotations: list of annotation file names
    :param outdir: str directory path in which to write revised annotations
    :param tmpdir: str path to pipeline's temporary directory
    :param seq_ident: int sequence identity threshold for matching duplicates
    :param seq_covg: int alignment coverage threshold for matching duplicates
    :return: list of paths to deduplicated annotation files
    """

    dedupe_tmp = os.path.join(tmpdir,'onegene')
    try:
        os.mkdir(dedupe_tmp)
    except:
        sys.exit("Could not create temporary directory in " + tmpdir)

    ref_cdss_all_fp = os.path.join(dedupe_tmp,'ref_cdss.faa')
    ref_cdss_unique_fp = os.path.join(dedupe_tmp,'ref_cdss_unique.faa')

    ann_sources = dict()
    ref_named_cds_count = {}
    with open(ref_cdss_all_fp,'w') as ref_cdss:
        for ref in annotations:
            ref_name = os.path.basename(os.path.splitext(ref)[0])
            ann_sources[ref_name] = ref
            ref_named_cds_count[ref_name] = extractor.fastaFromGbk(
                ref,
                out_cds = ref_cdss,
                out_genome = os.devnull,
                identify = lambda f: ':'.join([
                    ref_name,
                    extractor.get_ltag(f),
                    extractor.get_gene(f)
                    ])
            )
    if main_ref is None:
        main_ref = max(ref_named_cds_count, key=lambda _:ref_named_cds_count[_])

    rep_dict = CDHIT.run(
        nproc = 1,
        fasta = ref_cdss_all_fp,
        seq_ident = seq_ident,
        seq_covg = seq_covg,
        out = ref_cdss_unique_fp
    )


    # look for existing generic names in our gene names,
    # which are the third token in the colon-delimited string
    # Ref:locus_tag:gene_name
    existing_generigenes_fasta = os.path.join(dedupe_tmp,
                                              'preexisting_generigenes.fasta')
    extractor.subset_fasta(
        ref_cdss_all_fp,
        outseq = existing_generigenes_fasta,
        match = designator.is_unannotated,
        identify = lambda _: _.split(':')[2]
    )
    increment = designator.find_next_increment(existing_generigenes_fasta)


    subs = defaultdict(dict)
    subs_report = ''
    for rep in rep_dict.keys():
        subs, subs_report, increment = name_cluster(
            main_ref,
            [rep] + rep_dict[rep],
            increment,
            subs,
            subs_report
        )

    with open(os.path.join(outdir,'unifications.tsv'),'w') as report:
        report.write(subs_report)

    for ref in subs.keys():
        revised_records = []
        for record in SeqIO.parse(ann_sources[ref],'genbank'):
            revised = []
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers.keys() and \
                   feature.qualifiers['locus_tag'][0] in subs[ref].keys():
                    if 'gene' in feature.qualifiers.keys():
                        designator.append_qualifier(
                            feature.qualifiers,
                            'gene_synonym',
                            feature.qualifiers['gene'][0]
                        )
                        feature.qualifiers['gene'][0] = subs[ref][feature.qualifiers['locus_tag'][0]]
                    else:
                        feature.qualifiers['gene'] = [ subs[ref][feature.qualifiers['locus_tag'][0]] ]
                revised.append(feature)
            record.features = revised
            revised_records.append(record)
        ann_sources[ref] = os.path.abspath(os.path.join(outdir,ref + '.gbk'))
        SeqIO.write(revised_records, ann_sources[ref], format='genbank')

    return list(ann_sources.values())

def name_cluster(main_ref, cluster, increment, subs, subs_report):
    """
    choose a single gene name to apply to all members of a cluster.
    If there is one unique name and the rest of the cluster members are unnamed
    (that is, either an empty string or just the locus_tag), the existing name is
    used. If there are multiple unique names, or all members are unnamed, then a
    generic gene name is issued and incremented.

    :param main_ref: str name of reference for which to propagate locus tags as gene names
    :param cluster: list of strings representing cluster members,
                    in the form Ref:locus_tag:gene_name
    :param increment: last used number from the generic genes counter
    :param subs: nested dictionary giving the name to substitute for a locus tag in
                 a reference. Ex: subs['H37Rv']['Rv1164'] = 'ORF0001'
    :param subs_report: str in which to append tab-delimited summary of substitutions
                        in the following format:
                        Reference  locus_tag  original_gene_name  new_generic_name
    :returns: updated increment, subs, and subs_report
    """

    cluster_dict = defaultdict(list)
    cluster_list = []
    ref_gene_names = []
    for element in cluster:
        ref, ltag, gene_name = element.split(':')
        cluster_dict[ref].append((ltag, gene_name))
        cluster_list.append((ref, ltag, gene_name))
        if gene_name and gene_name != ltag:
            ref_gene_names.append((ref, gene_name))


    #
    # only one unique name exists in the cluster
    # => use it for all genes in all refs
    #
    if len(ref_gene_names) >= 1 and all([gene==ref_gene_names[0][1] for ref, gene in ref_gene_names]):
        # Make sure there are at least some unnamed genes here
        if len(ref_gene_names) != len(cluster):
            name = ref_gene_names[0][1]
        # Nothing to do; all genes uniformly named
        else:
            cluster_list = []
    #
    # no names exist
    # => use the primary reference's locus tag as the name for all genes
    #    (so long as the primary reference only has a single copy; otherwise, assign generic)
    #
    elif not ref_gene_names:
        # Make sure there is more than one gene in this case
        if len(cluster) > 1:
            primaries = [(ref, gene) for ref, ltag, gene in cluster_list if ref==main_ref]
            if len(primaries) == 1:
                name = primaries[0][1]
            else:
                (name, increment) = designator.assign_orf_id(increment)
        # We won't assign a generic name for just a single gene.
        else:
            cluster_list = []
    #
    # no names exist or multiple unique names exist, but they all come from the same reference annotation
    # => assign a generic name and use it for all genes in all refs
    #
    elif all([ref==ref_gene_names[0][0] for ref, gene in ref_gene_names]):
        # Make sure there is more than one gene in this case
        if len(cluster) > 1:
            (name, increment) = designator.assign_orf_id(increment)
        # We won't assign a generic name for just a single gene.
        else:
            cluster_list = []
    #
    # multiple unique names exist across references
    # => process individually
    #
    else:
        for ref in cluster_dict:
            subs, subs_report, increment = name_cluster(
                main_ref,
                [_ for _ in cluster if _.startswith(f"{ref}:")],
                increment,
                subs,
                subs_report,
            )
        cluster_list = []


    for ref, ltag, og_gene_name in cluster_list:
        subs[ref][ltag] = name
        subs_report += '\t'.join([
            ref,
            ltag,
            og_gene_name,
            name]) + '\n'

    return subs, subs_report, increment
