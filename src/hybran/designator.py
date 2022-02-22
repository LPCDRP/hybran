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

def find_next_increment(fasta):
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
        return int(last_orf.replace(generic_orf_prefix, '')) + 1
    else:
        return 1

#
# These can be applied to sequence record IDs to match the designated property
#

def is_unannotated(name):
    return name.startswith(generic_orf_prefix)

def is_reference(name):
    return not name.startswith((generic_orf_prefix,'L_','L2_'))
