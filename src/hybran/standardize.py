import os
import re
import warnings

from . import designator
from . import fileManager
from .bio import SeqIO


def main(args):
    annotations = fileManager.file_list(
        args.annotations,
        file_type='genbank'
    )

    designator.generic_orf_prefix[0]=args.orf_prefix
    designator.ref_orf_prefix[0] = f"REF{args.orf_prefix}X"

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit("Could not create output directory " + args.output)

    if not os.path.isfile(args.unifications_file):
        sys.exit(f"Invalid unifications file {args.unifications_file}")
    generics = load_reference_generics(args.unifications_file)

    for annotation in annotations:
        outfile = os.path.join(
            args.output,
            os.path.basename(annotation),
        )
        records = list(SeqIO.parse(annotation, 'genbank'))
        for record in records:
            for feature in record.features:
                standardize(feature, generics)
        with open(outfile, 'w') as outfile_handle:
            SeqIO.write(records, outfile_handle, 'genbank')


def standardize(feature, generics):
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
        if fusion and (
                designator.has_unannotated_component(gene)
                or designator.has_uniref_component(gene)
        ):
            components = gene.split('::')
            for i in range(len(components)):
                if components[i] in generics:
                    components[i] = generics[components[i]]
                    if 'inference' in feature.qualifiers:
                        update_inferences(
                            feature.qualifiers['inference'],
                            generics,
                        )
                elif designator.is_unannotated(components[i]):
                    # HGNC nomenclature for fusion with an unknown gene
                    components[i] = '?'
            feature.qualifiers['gene'][0] = '::'.join(components)
        elif not fusion and (
                designator.is_unannotated(feature.qualifiers['gene'][0])
                or designator.is_uniref(feature.qualifiers['gene'][0])
        ):
            if 'gene_synonym' in feature.qualifiers:
                if 'inference' in feature.qualifiers:
                    update_inferences(
                        feature.qualifiers['inference'],
                        {
                            feature.qualifiers['gene'][0]:feature.qualifiers['gene_synonym'][-1]
                        },
                    )
                feature.qualifiers['gene'][0] = feature.qualifiers['gene_synonym'].pop()
                if len(feature.qualifiers['gene_synonym']) == 0:
                    del feature.qualifiers['gene_synonym']
            elif gene in generics:
                feature.qualifiers['gene'][0] = generics[gene]
                if 'inference' in feature.qualifiers:
                    update_inferences(
                        feature.qualifiers['inference'],
                        generics,
                    )
            else:
                del feature.qualifiers['gene']


def load_reference_generics(infile):
    """
    Load the unified-refs/unifications.tsv file into a dictionary.
    """
    generics = {}
    with open(infile, 'r') as unitable:
        for line in unitable:
            [ref, ltag, gene, generic_name] = line.strip().split('\t')
            generics[generic_name] = gene
    return generics

def update_inferences(inference_list, generics):
    for i in range(len(inference_list)):
        try:
            inference_list[i] = re.sub(
                designator.ref_orf_prefix[0] + r"\d+",
                lambda _: generics[_[0]],
                inference_list[i],
            )
        except KeyError:
            warnings.warn(f"Did not find matching generic name to replace that in the inference: {inference_list[i]}")
