from io import StringIO

from Bio import SeqIO


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

    with StringIO() as fasta:
        # Check for existing instances of this locus tag prefix
        SeqIO.convert(gbk, "genbank", fasta, "fasta")
        increment = find_next_increment(fasta, prefix=prefix)

    old_to_new = dict()
    for record in SeqIO.parse(gbk, "genbank"):
        output_records.append(record)
        for feature in record.features:
            if 'locus_tag' not in feature.qualifiers.keys():
                continue
            elif feature.qualifiers['locus_tag'][0].startswith(prefix+'_'):
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
                    new_locus_tag = '_'.join([prefix, "%04g" % (increment)])
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

def find_next_increment(fasta, prefix=generic_orf_prefix):
    """
    Based on a given unannotated FASTA, identifies the highest
    numbered ORF increment in a list of numbered features and returns
    the next number to use.

    :param fasta: FASTA file name
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
