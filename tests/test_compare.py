from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pytest

from hybran import compare


@pytest.mark.parametrize('case', [
    'two_word_code',
    'multiple_codes',
    'not_relevant',
    'ambiguous_residues',
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
        'ambiguous_residues': (
            "frameshifted; too many ambiguous residues; Derived "
            "by automated computational analysis using gene prediction "
            "method: Protein Homology."
        ),
    }
    expected = {
        'two_word_code': 'internal_stop',
        'multiple_codes': 'internal_stop;incomplete',
        'not_relevant': '.',
        'ambiguous_residues': 'frameshifted;too_many_ambiguous_residues',
    }

    assert compare.pgap_np(
        SeqFeature(qualifiers={'note': [ note[case] ]})
    ) == expected[case]


@pytest.mark.parametrize('case', [
    'colo_conflict_fusion',
    'colo_conflict_fusion_swapped',
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
    swapped = False
    if case.endswith('_swapped'):
        swapped = True
        case = case.replace('_swapped','')

    if not swapped:
        f1_list, f2_list = inputs[case]
    else:
        f2_list, f1_list = inputs[case]

    expected = {
        'colo_conflict_fusion': (
            # colocated
            [
                ( inputs[case][0][0], inputs[case][1][0] ),
                ( inputs[case][0][1], inputs[case][1][2] ),
                ( inputs[case][0][2], inputs[case][1][3] ),
            ],
            # conflicting
            [],
            # unique in 1
            [],
            # unique in 2
            [
                inputs[case][1][1],
            ]
        ),
    }

    # comparison results should be identical regardless of which is called genome 1
    if swapped:
        expected[case] = list(expected[case])
        expected[case][0] = [(y, x) for x,y in expected[case][0]]
        expected[case][1] = [(y, x) for x,y in expected[case][1]]
        expected[case][2], expected[case][3] = expected[case][3], expected[case][2]
        expected[case] = tuple(expected[case])

    assert compare.compare(f1_list, f2_list)[0:4] == expected[case]
