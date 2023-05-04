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

    dedupe(
        annotations=fileManager.file_list(args.annotations, file_type="genbank"),
        outdir=args.output,
        tmpdir=config.hybran_tmp_dir,
        seq_ident=args.identity_threshold,
        seq_covg=args.coverage_threshold,
    )

def dedupe(annotations, outdir, tmpdir, seq_ident=99, seq_covg=99):
    """"
    Identify duplicate genes among the reference annotation and assign all instances a
    single name. If dupicates had different actual gene names, these original names are
    stored in the gene_synonym field.

    :param annotations: list of annotation file names
    :param outdir: str directory path in which to write revised annotations
    :param tmpdir: str path to pipeline's temporary directory
    :param seq_ident: int sequence identity threshold for matching duplicates
    :param seq_covg: int alignment coverage threshold for matching duplicates
    :return: list of paths to deduplicated annotation files
    """

    dedupe_tmp = os.path.join(tmpdir,'ref-dedupe')
    try:
        os.mkdir(dedupe_tmp)
    except:
        sys.exit("Could not create temporary directory in " + tmpdir)

    ref_cdss_all_fp = os.path.join(dedupe_tmp,'ref_cdss.faa')
    ref_cdss_unique_fp = os.path.join(dedupe_tmp,'ref_cdss_unique.faa')

    ann_sources = dict()
    with open(ref_cdss_all_fp,'w') as ref_cdss:
        for ref in annotations:
            ref_name = os.path.basename(os.path.splitext(ref)[0])
            ann_sources[ref_name] = ref
            extractor.fastaFromGbk(
                ref,
                out_cds = ref_cdss,
                out_genome = os.devnull,
                identify = lambda f: ':'.join([
                    ref_name,
                    extractor.get_ltag(f),
                    extractor.get_gene(f)
                    ])
            )
    
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
            [rep] + rep_dict[rep],
            increment,
            subs,
            subs_report
        )

    with open(os.path.join(outdir,'duplicates.tsv'),'w') as report:
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

def name_cluster(cluster, increment, subs, subs_report):
    """
    choose a single gene name to apply to all members of a cluster.
    If there is one unique name and the rest of the cluster members are unnamed
    (that is, either an empty string or just the locus_tag), the existing name is
    used. If there are multiple unique names, or all members are unnamed, then a
    generic gene name is issued and incremented.

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
    for element in cluster:
        ref, ltag, gene_name = element.split(':')
        cluster_dict[ref].append((ltag, gene_name))

    for ref in cluster_dict:
        if len(cluster_dict[ref]) > 1:
            # see how many actual names we have
            actual_gene_names = [gene for ltag, gene in cluster_dict[ref] if gene and gene != ltag]

            if len(set(actual_gene_names)) == 1:
                name = actual_gene_names[0]
            else:
                (name, increment) = designator.assign_orf_id(increment)

            for ltag, gene in cluster_dict[ref]:
                subs[ref][ltag] = name
                subs_report += '\t'.join(
                    [ref,
                     ltag,
                     gene,
                     name]) + '\n'

    return subs, subs_report, increment
