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

    for annotation in annotations:
        outfile = os.path.join(
            args.output,
            os.path.basename(annotation),
        )
        records = list(SeqIO.parse(annotation, 'genbank'))
        for record in records:
            for feature in record.features:
                standardize(feature)
        with open(outfile, 'w') as outfile_handle:
            SeqIO.write(records, outfile_handle, 'genbank')


def standardize(feature):
    """
    remove generic gene names.
    If a gene_synonym is available (as from an ab initio call or
    a RATT-transferred annotation of a reference paralog), use it.
    If not, remove the `gene` qualifier altogether.

    :param feature: SeqFeature object
    """
    if ('gene' in feature.qualifiers
        and designator.is_unannotated(feature.qualifiers['gene'][0])
        ):
        if 'gene_synonym' in feature.qualifiers:
            feature.qualifiers['gene'][0] = feature.qualifiers['gene_synonym'].pop()
        else:
            del feature.qualifiers['gene']
