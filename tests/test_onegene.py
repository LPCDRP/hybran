from collections import defaultdict

import pytest

from hybran import designator, onegene

@pytest.mark.parametrize("situation", [
    'singleref_multiple_nonames',
    'singleref_one_named',
    'multiref_single_orthologs',
    'multiref_mixed_bag',
])
def test_name_cluster(situation):

    inputs = {
        'singleref_multiple_nonames': [
            'H37Rv:Rv1369c:Rv1369c',
            'H37Rv:Rv1756c:Rv1756c',
            'H37Rv:Rv1764:Rv1764',
        ],
        'singleref_one_named': [
            'H37Rv:Rv1454c:qor',
        ],
        'multiref_single_orthologs': [
            'Ref1:R1_01234:R1_01234',
            'Ref2:R2_01234:R2_01234',
        ],
        'multiref_mixed_bag': [
            'Ref1:R1_01234:R1_01234',
            'Ref1:R1_05678:',
            'Ref2:R2_01234:R2_01234',
            'Ref3:R3_01234:genA',
            'Ref3:R3_06789:R3_06789',
        ],
    }

    designator.generic_orf_prefix = ['ORF']
    expected_subs = defaultdict(dict)
    expected = {
        'singleref_multiple_nonames': [
            {'H37Rv': {
                'Rv1369c':'ORF0001',
                'Rv1756c':'ORF0001',
                'Rv1764':'ORF0001',
            }},
            ("H37Rv\tRv1369c\tRv1369c\tORF0001\n"
             "H37Rv\tRv1756c\tRv1756c\tORF0001\n"
             "H37Rv\tRv1764\tRv1764\tORF0001\n"),
            2,
        ],
        'singleref_one_named': [
            {},
            "",
            1,
        ],
        'multiref_single_orthologs': [
            {},
            "",
            1,
        ],
        'multiref_mixed_bag': [
            {'Ref1': {
                'R1_01234':'ORF0001',
                'R1_05678':'ORF0001',
            },
             'Ref3': {
                 'R3_01234':'genA',
                 'R3_06789':'genA',
             },
             },
            ("Ref1\tR1_01234\tR1_01234\tORF0001\n"
             "Ref1\tR1_05678\t\tORF0001\n"
             "Ref3\tR3_01234\tgenA\tgenA\n"
             "Ref3\tR3_06789\tR3_06789\tgenA\n"
             ),
            2,
        ],
    }
    expected_subs.update(expected[situation][0])
    expected[situation][0] = expected_subs
    assert onegene.name_cluster(
        cluster=inputs[situation],
        increment=1,
        subs=defaultdict(dict),
        subs_report='',
    ) == tuple(expected[situation])
