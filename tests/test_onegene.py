from collections import defaultdict

import pytest

from hybran import designator, onegene

@pytest.mark.parametrize("situation", [
    'singleref_multiple_nonames',
    'singleref_one_named',
    'multiref_single_orthologs',
    'multiref_mixed_bag',
    'multiref_orthologous_named_duplicates',
])
def test_name_cluster(situation):

    inputs = {
        'singleref_multiple_nonames': [
            'Ref1:Rv1369c:Rv1369c',
            'Ref1:Rv1756c:Rv1756c',
            'Ref1:Rv1764:Rv1764',
        ],
        'singleref_one_named': [
            'Ref1:Rv1454c:qor',
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
        'multiref_orthologous_named_duplicates': [
            'Ref1:R1_01234:genA',
            'Ref1:R1_01235:genB',
            'Ref2:R2_01234:genA',
            'Ref2:R2_01235:genB',
        ],
    }

    designator.generic_orf_prefix = ['ORF']
    expected_subs = defaultdict(dict)
    expected = {
        'singleref_multiple_nonames': [
            {'Ref1': {
                'Rv1369c':'ORF0001',
                'Rv1756c':'ORF0001',
                'Rv1764':'ORF0001',
            }},
            ("Ref1\tRv1369c\tRv1369c\tORF0001\n"
             "Ref1\tRv1756c\tRv1756c\tORF0001\n"
             "Ref1\tRv1764\tRv1764\tORF0001\n"),
            2,
        ],
        'singleref_one_named': [
            {},
            "",
            1,
        ],
        'multiref_single_orthologs': [
            {'Ref1': {
                'R1_01234':'R1_01234',
            },
             'Ref2': {
                'R2_01234':'R1_01234',
            },
            },
            ("Ref1\tR1_01234\tR1_01234\tR1_01234\n"
             "Ref2\tR2_01234\tR2_01234\tR1_01234\n"),
            1,
        ],
        'multiref_mixed_bag': [
            {'Ref1': {
                'R1_01234':'genA',
                'R1_05678':'genA',
            },
             'Ref2': {
                 'R2_01234': 'genA',
             },
             'Ref3': {
                 'R3_01234':'genA',
                 'R3_06789':'genA',
             },
             },
            ("Ref1\tR1_01234\tR1_01234\tgenA\n"
             "Ref1\tR1_05678\t\tgenA\n"
             "Ref2\tR2_01234\tR2_01234\tgenA\n"
             "Ref3\tR3_01234\tgenA\tgenA\n"
             "Ref3\tR3_06789\tR3_06789\tgenA\n"
             ),
            1,
        ],
        'multiref_orthologous_named_duplicates': [
            {'Ref1': {
                'R1_01234':'ORF0001',
                'R1_01235':'ORF0001',
            },
             'Ref2': {
                'R2_01234':'ORF0001',
                'R2_01235':'ORF0001',
             },
             },
            ("Ref1\tR1_01234\tgenA\tORF0001\n"
             "Ref1\tR1_01235\tgenB\tORF0001\n"
             "Ref2\tR2_01234\tgenA\tORF0001\n"
             "Ref2\tR2_01235\tgenB\tORF0001\n"
             ),
            2,
        ],
    }
    expected_subs.update(expected[situation][0])
    expected[situation][0] = expected_subs
    assert onegene.name_cluster(
        main_ref="Ref1",
        cluster=inputs[situation],
        increment=1,
        subs=defaultdict(dict),
        subs_report='',
    ) == tuple(expected[situation])
