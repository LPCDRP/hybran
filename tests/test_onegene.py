from collections import defaultdict

import pytest

from hybran import designator, onegene

@pytest.mark.parametrize("situation", [
    'singleref_multiple_nonames',
    'singleref_one_named',
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
    }
    expected_subs.update(expected[situation][0])
    expected[situation][0] = expected_subs
    assert onegene.name_cluster(
        cluster=inputs[situation],
        increment=1,
        subs=defaultdict(dict),
        subs_report='',
    ) == tuple(expected[situation])
