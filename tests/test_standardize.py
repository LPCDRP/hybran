from copy import deepcopy

from Bio.SeqFeature import SeqFeature
import pytest

from hybran import designator
from hybran import standardize

@pytest.mark.parametrize('gene', [
    'named',
    'unnamed',
    'unnamed_with_syn',
    'refonly_unnamed_with_syn',
    'named_unnamed_with_syn',
    'named_unnamed',
    'named_unlisted',
    'unnamed_named',
    'named_unnamed_named',
    'unnamed_named,unnamed',
])
def test_standardize(gene):
    inputs = {
        'named': SeqFeature(
            qualifiers={'gene':["MBT_01234"]}
        ),
        'unnamed': SeqFeature(
            qualifiers={'gene':["ORF2345"]}
        ),
        'unnamed_with_syn': SeqFeature(
            qualifiers={'gene':["ORF2395"], 'gene_synonym':['MBT_54632']}
        ),
        'refonly_unnamed_with_syn': SeqFeature(
            qualifiers={'gene':["ORF1111"], 'gene_synonym':['MBT_54632']}
        ),
        'named_unnamed': SeqFeature(
            qualifiers={'gene':["MBT_01234::ORF2346"]}
        ),
        'named_unnamed_with_syn': SeqFeature(
            qualifiers={'gene':["MBT_01234::ORF4356"], 'gene_synonym':['mgv47']}
        ),
        'named_unlisted': SeqFeature(
            qualifiers={'gene':["MBT_01234::ORF0000"]}
        ),
        'unnamed_named': SeqFeature(
            qualifiers={'gene':["ORF4356::MBT_01234"]}
        ),
        'named_unnamed_named': SeqFeature(
            qualifiers={'gene':["MBT_01234::ORF4545::MBT_4356"]}
        ),
        'unnamed_named,unnamed': SeqFeature(
            qualifiers={'gene':["ORF4356::MBT_01234::ORF5467"]}
        ),
    }
    expected = {
        'named': SeqFeature(
            qualifiers={'gene':["MBT_01234"]}
        ),
        'unnamed': SeqFeature(
            qualifiers={'gene':["MBT_4356"]}
        ),
        'unnamed_with_syn': SeqFeature(
            qualifiers={'gene':["MBT_54632"]},
        ),
        'refonly_unnamed_with_syn': SeqFeature(
            qualifiers={'gene_synonym':['MBT_54632']}
        ),
        'named_unnamed': SeqFeature(
            qualifiers={'gene':["MBT_01234::MBT_3456"]}
        ),
        'named_unnamed_with_syn': SeqFeature(
            qualifiers={'gene':["MBT_01234::MBT_5678"], 'gene_synonym':['mgv47']}
        ),
        'named_unlisted': SeqFeature(
            qualifiers={'gene':["MBT_01234::?"]}
        ),
        'unnamed_named': SeqFeature(
            qualifiers={'gene':["MBT_5678::MBT_01234"]}
        ),
        'named_unnamed_named': SeqFeature(
            qualifiers={'gene':["MBT_01234::MBT_9876::MBT_4356"]}
        ),
        'unnamed_named,unnamed': SeqFeature(
            qualifiers={'gene':["MBT_5678::MBT_01234::MBT_4359"]}
        ),
    }
    generics = {
        'ORF2345': 'MBT_4356',
        'ORF2495': 'MBT_5467',
        'ORF2395': 'MBT_5457',
        'ORF2346': 'MBT_3456',
        'ORF4356': 'MBT_5678',
        'ORF4545': 'MBT_9876',
        'ORF5467': 'MBT_4359',
    }
    designator.generic_orf_prefix[0]= 'ORF'
    ref_only = gene.startswith('refonly')

    standardize.standardize(inputs[gene], generics=generics, ref_only=ref_only)
    assert inputs[gene].qualifiers == expected[gene].qualifiers
