from collections import defaultdict
import re

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from . import extractor

generic_orf_prefix = 'ORF'


def append_qualifier(qualifiers, qual_name, qual_value):
    """
    Append to a list or create a new one if it doesn't exist

    :param qualifiers: dict (pass feature.qualifiers here)
    :param qual_name: str feature qualifier key
    :param qual_value: str
    :return: None (original dictionary is modified)
    """
    if qual_name in qualifiers.keys():
        qualifiers[qual_name].append(qual_value)
    else:
        qualifiers[qual_name] = [qual_value]

def assign_locus_tags(gbk, prefix):

    output_records = []

    # Check for existing instances of this locus tag prefix
    # TODO - refactoring opportunity with find_next_increment()
    ltags = []
    n_digits = 5
    delim = '_'
    for ltag in extractor.gene_dict(gbk, cds_only=False).keys():
        if ltag.startswith(prefix):
            ltags.append(ltag)
    if ltags:
        # some annotations don't use the underscore to separate the prefix,
        if not ltags[0].startswith(prefix + '_'):
            delim = ''
        ltag_numbers = [re.sub(rf'^{prefix}{delim}(\d+).*',r'\1', _) for _ in ltags]
        n_digits = len(ltag_numbers[0])
        increment = int(sorted(ltag_numbers)[-1]) + 1
    else:
        increment = 1

    old_to_new = dict()
    for record in SeqIO.parse(gbk, "genbank"):
        output_records.append(record)
        for feature in record.features:
            if 'locus_tag' not in feature.qualifiers.keys():
                continue
            elif feature.qualifiers['locus_tag'][0].startswith(prefix+delim):
                continue
            else:
                orig_locus_tag = feature.qualifiers['locus_tag'][0]
                # For the case of multiple features that must have the
                # same locus tag (like gene + mRNA + CDS + 5'UTR + ...)
                # We want to use the same new locus tag for all of these.
                if orig_locus_tag in old_to_new.keys():
                    feature.qualifiers['locus_tag'][0] = \
                        old_to_new[orig_locus_tag]
                else:
                    new_locus_tag = delim.join([prefix, f"%0{n_digits}g" % (increment)])
                    increment += 1
                    old_to_new[orig_locus_tag] = new_locus_tag
                    feature.qualifiers['locus_tag'][0] = new_locus_tag

    SeqIO.write(output_records, gbk, "genbank")

def assign_orf_id(increment):
    """
    Format an ID string for a feature at the given count.

    :param increment: int count
    :returns:
       - str ID string
       - incremented counter
    """
    orf_id = generic_orf_prefix + "%04g" % (increment)
    increment += 1
    return orf_id, increment

def create_gene_entries(gbk):

    output_records = []
    gene_positions = defaultdict(list)
    genes = dict()
    for record in SeqIO.parse(gbk, "genbank"):
        updated_record_features = []
        for f in record.features:
            if f.type == 'gene':
                genes[f.qualifiers['locus_tag'][0]] = f
            elif 'locus_tag' in f.qualifiers.keys():
                ltag = f.qualifiers['locus_tag'][0]
                if ltag not in genes.keys():
                    genes[ltag] = SeqFeature(
                        FeatureLocation(
                            f.location.start,
                            f.location.end,
                            f.location.strand,
                        ),
                        type = 'gene',
                        qualifiers = dict(
                            locus_tag=ltag,
                            gene=extractor.get_gene(f)
                        )
                    )
                    updated_record_features.append(genes[ltag])
                else:
                    genes[ltag].location = FeatureLocation(
                        genes[ltag].location.start,
                        f.location.end,
                        genes[ltag].location.strand
                    )
                if 'pseudo' in f.qualifiers.keys():
                    genes[ltag].qualifiers['pseudo'] = ['']
            updated_record_features.append(f)
        record.features = updated_record_features
        output_records.append(record)

    SeqIO.write(output_records, gbk, format="genbank")

def find_next_increment(fasta, prefix=generic_orf_prefix):
    """
    Based on a given unannotated FASTA, identifies the highest
    numbered ORF increment in a list of numbered features and returns
    the next number to use.

    :param fasta: FASTA file name
    :param format: input file format
    :return: int
    """
    records = []
    for record in SeqIO.parse(fasta, 'fasta'):
        records.append(record.id)
    if records:
        last_orf = sorted(records)[-1]
        return int(last_orf.replace(prefix, '')) + 1
    else:
        return 1

#
# These can be applied to sequence record IDs to match the designated property
#

def is_unannotated(name):
    return name.startswith(generic_orf_prefix)

def is_reference(name):
    return not name.startswith((generic_orf_prefix,'L_','L2_'))

def is_raw_ltag(name):
    return name.startswith(('L_','L2_'))
