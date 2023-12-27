from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pytest

from hybran import compare


@pytest.mark.parametrize('case', [
    'two_word_code',
    'multiple_codes',
    'not_relevant',
])
def test_pgap_np(case):
    note = {
        'two_word_code': (
            "internal stop; Derived by automated computational "
            "analysis using gene prediction method: Protein Homology. "
            "GO_function: GO:0003677 - DNA binding [Evidence IEA]; "
            "GO_function: GO:0004519 - endonuclease activity [Evidence "
            "IEA]"
        ),
        'multiple_codes': (
            "internal stop; incomplete; partial on complete "
            "genome; missing C-terminus; Derived by automated "
            "computational analysis using gene prediction method: "
            "Protein Homology."
        ),
        'not_relevant': (
            "Essential for efficient processing of 16S rRNA; "
            "Derived by automated computational analysis using gene "
            "prediction method: Protein Homology. "
            "GO_function: GO:0003723 - RNA binding [Evidence IEA]; "
            "GO_process: GO:0006364 - rRNA processing [Evidence IEA]"
        ),
    }
    expected = {
        'two_word_code': 'internal_stop',
        'multiple_codes': 'internal_stop;incomplete',
        'not_relevant': '.',
    }

    assert compare.pgap_np(
        SeqFeature(qualifiers={'note': [ note[case] ]})
    ) == expected[case]


@pytest.mark.parametrize('case', [
    'colo_conflict_fusion',
])
def test_compare(case):
    inputs = {
        'colo_conflict_fusion': [
            [
                # pgap
                SeqFeature(
                    FeatureLocation(959902, 960772, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_000905'], 'gene':['trxA']}
                ),
                SeqFeature(
                    FeatureLocation(960843, 961506, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_000906']}
                ),
                SeqFeature(
                    FeatureLocation(961502, 961973, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_000907']}
                ),
            ], [
                # hybran
                SeqFeature(
                    FeatureLocation(959902, 960772, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_00926'], 'gene':['PA4061']}
                ),
                SeqFeature(
                    FeatureLocation(960843, 961122, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_00927'], 'gene':['PA4060']}
                ),
                SeqFeature(
                    FeatureLocation(960843, 961506, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_00928'], 'gene':['PA4059::PA4060']}
                ),
                SeqFeature(
                    FeatureLocation(961502, 961973, -1), type='CDS',
                    qualifiers={'locus_tag':['PAK_00929'], 'gene':['PA4058']}
                ),
            ]
        ],
    }
    f1_list, f2_list = inputs[case]

    expected = {
        'colo_conflict_fusion': (
            # colocated
            [
                [
                    'PAK_000905', 'trxA', 0, '.',
                    'PAK_00926', 'PA4061', 0, '.',
                    ExactPosition(959902), ExactPosition(960772), -1
                ], [
                    'PAK_000906', 'PAK_000906', 0, '.',
                    'PAK_00928', 'PA4059::PA4060', 0, '.',
                    ExactPosition(960843), ExactPosition(961506),  -1
                ], [
                    'PAK_000907', 'PAK_000907', 0, '.',
                    'PAK_00929', 'PA4058', 0, '.',
                    ExactPosition(961502), ExactPosition(961973), -1
                ],
            ],
            # conflicting
            [],
            # unique
            []
        ),
    }

    assert compare.compare(f1_list, f2_list) == expected[case]
