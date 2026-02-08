from Bio.SeqFeature import SeqFeature, SimpleLocation
import pytest

from hybran import defuse


@pytest.mark.parametrize('situation', [
    'whole_gene_fusion', # https://gitlab.com/LPCDRP/hybran/-/issues/81
])
def test_defuse(situation):
    if situation == 'whole_gene_fusion':
        feature_list = [
            SeqFeature(
                SimpleLocation(2360919, 2362219, -1),
                type='CDS',
                qualifiers={'locus_tag':['1-0006_02304'], 'gene':['PE_PGRS36'], 'pseudo':['']}
            ),
            SeqFeature(
                SimpleLocation(2360919, 2362392, -1),
                type='CDS',
                qualifiers={'locus_tag':['1-0006_02305'], 'gene':['PE21::PE_PGRS36']}
            )
        ]
        fusions_by_ltag = {
            '1-0006_02305': {
                'type': 'whole',
                'gene': 'PE21::PE_PGRS36',
                'strand': '-1',
                'start': '2360920',
                'end': '2362392',
                'components': [
                    SeqFeature(
                        SimpleLocation(2362218, 2362392, strand=-1),
                        type='CDS',
                        qualifiers={'gene':['PE21'], 'pseudo':['']}
                    ),
                    SeqFeature(
                        SimpleLocation(2360919, 2362219, strand=-1),
                        type='CDS',
                        qualifiers={'gene':['PE_PGRS36'], 'locus_tag':['1-0006_02304']}
                    )
                ],
            }
        }
        expected = [
            SeqFeature(
                SimpleLocation(2360919, 2362219, -1),
                type='CDS',
                qualifiers={'locus_tag':['1-0006_02304'], 'gene':['PE_PGRS36'], 'pseudo':['']}
            ),
            SeqFeature(
                SimpleLocation(2362218, 2362392, -1),
                type='CDS',
                qualifiers={'locus_tag':['1-0006_02305'], 'gene':['PE21'], 'pseudo':['']}
            )
        ]

    assert defuse.defuse(
        feature_list,
        fusions_by_ltag,
    ) == expected
