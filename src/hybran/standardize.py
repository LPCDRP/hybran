import os

from Bio import SeqIO

from . import designator
from . import fileManager


def main(args):
    annotations = fileManager.file_list(
        args.annotations,
        file_type='genbank'
    )

    designator.generic_orf_prefix[0]=args.orf_prefix

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit("Could not create output directory " + args.output)

    if not os.path.isfile(args.duplicates_file):
        sys.exit(f"Invalid duplicates file {args.duplicates_file}")
    duplicates = load_reference_duplicates(args.duplicates_file)

    for annotation in annotations:
        outfile = os.path.join(
            args.output,
            os.path.basename(annotation),
        )
        records = list(SeqIO.parse(annotation, 'genbank'))
        for record in records:
            for feature in record.features:
                standardize(feature, duplicates)
        with open(outfile, 'w') as outfile_handle:
            SeqIO.write(records, outfile_handle, 'genbank')


def standardize(feature, duplicates):
    """
    remove generic gene names.
    If a gene_synonym is available (as from an ab initio call or
    a RATT-transferred annotation of a reference paralog), use it.
    If not, remove the `gene` qualifier altogether.

    :param feature: SeqFeature object
    """
    if 'gene' in feature.qualifiers:
        gene = feature.qualifiers['gene'][0]
        fusion = '::' in gene
        if fusion and designator.has_unannotated_component(gene):
            components = gene.split('::')
            for i in range(len(components)):
                if designator.is_unannotated(components[i]):
                    if components[i] in duplicates:
                        components[i] = duplicates[components[i]]
                    else:
                        # HGNC nomenclature for fusion with an unknown gene
                        components[i] = '?'
            feature.qualifiers['gene'][0] = '::'.join(components)
        elif not fusion and designator.is_unannotated(feature.qualifiers['gene'][0]):
            if 'gene_synonym' in feature.qualifiers:
                feature.qualifiers['gene'][0] = feature.qualifiers['gene_synonym'].pop()
                if len(feature.qualifiers['gene_synonym']) == 0:
                    del feature.qualifiers['gene_synonym']
            elif gene in duplicates:
                feature.qualifiers['gene'][0] = duplicates[gene]
            else:
                del feature.qualifiers['gene']


def load_reference_duplicates(infile):
    """
    Load the deduped-refs/duplicates.tsv file into a dictionary.
    """
    duplicates = {}
    with open(infile, 'r') as dupe_table:
        for line in dupe_table:
            [ref, ltag, gene, generic_name] = line.strip().split('\t')
            duplicates[generic_name] = gene
    return duplicates
