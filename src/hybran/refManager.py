from collections import defaultdict
import os

from Bio import SeqIO

from . import CDHIT, extractor, parseClustering


def dedupe(annotations, outdir, seq_ident=99, seq_covg=99):
    """"
    Identify duplicate genes among the reference annotation and assign all instances a
    single name. If dupicates had different actual gene names, these original names are
    stored in the gene_synonym field.

    :param annotations: list of annotation file names
    :param outdir: str directory path in which to write revised annotations
    :param seq_ident: int sequence identity threshold for matching duplicates
    :param seq_covg: int alignment coverage threshold for matching duplicates
    :return: list of paths to deduplicated annotation files
    """

    ref_cdss_all_fp = os.path.join(outdir,'ref_cdss.faa')
    ref_cdss_unique_fp = os.path.join(outdir,'ref_cdss_unique.faa')

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


    # TODO - this function currently doesn't work with our input since our seq id headers are not
    #        just the gene names. they're Ref:locus_tag:gene_name
    #increment = parseClustering.find_largest_mtb_increment(ref_cdss_all_fp)
    increment = 0

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
        print(subs_report, file=report)

    for ref in subs.keys():
        revised_records = []
        for record in SeqIO.parse(ann_sources[ref],'genbank'):
            revised = []
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers.keys() and \
                   feature.qualifiers['locus_tag'][0] in subs[ref].keys():
                    if feature.qualifiers['gene'][0]:
                        if 'gene_synonym' in feature.qualifiers.keys():
                            feature.qualifiers['gene_synonym'].append(
                                feature.qualifiers['gene'][0]
                            )
                        else:
                            feature.qualifiers['gene_synonym'] = feature.qualifiers['gene'].copy()
                    feature.qualifiers['gene'][0] = subs[ref][feature.qualifiers['locus_tag'][0]]
                revised.append(feature)
            record.features = revised
            revised_records.append(record)
        ann_sources[ref] = os.path.abspath(os.path.join(outdir,ref + '.gbk'))
        SeqIO.write(revised_records, ann_sources[ref], format='genbank')

    return list(ann_sources.values())

def name_cluster(cluster, increment, subs, subs_report):
    ref_names, ltags, gene_names = list(zip(*(_.split(':') for _ in cluster)))

    if len(cluster) == 1:
        name = gene_names[0]
    else:
        # see how many actual names we have
        actual_gene_names = [gene_names[i] for i in range(len(gene_names)) if gene_names[i] != ltags[i]]
        unnamed_genes = [gene_names[i] for i in range(len(gene_names)) if gene_names[i] == ltags[i]]

        if len(set(actual_gene_names)) == 1:
            name = actual_gene_names[0]
        else:
            name = 'MTB' + "%04g" % (1 + increment)
            increment += 1

        for i in range(len(cluster)):
            if gene_names[i] != name:
                subs[ref_names[i]][ltags[i]] = name
                subs_report += '\t'.join(
                    [ref_names[i],
                     ltags[i],
                     gene_names[i],
                     name]) + '\n'

    return subs, subs_report, increment
