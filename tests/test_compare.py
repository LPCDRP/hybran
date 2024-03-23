from Bio.SeqFeature import (
    SeqFeature,
    SimpleLocation,
    CompoundLocation,
    FeatureLocation,
    ExactPosition,
)
import pytest

from hybran import compare

from .data_features import *


@pytest.mark.parametrize('pair', [
    'two_pseudos_same_start',
    'same_stop',
    'inframe_pseudo_same_start_different_stop',
    'inframe_pseudo_different_start_different_stop',
    'inframe_pseudo_on_nonoverlapping_part',
    'overlapping_out_of_frame',
    'different_strand',
    'non_overlapping',
    'one_bp_apart',
    'compound_interval',
])
def test_overlap_inframe(pair):
    pairs = {
        'two_pseudos_same_start': (
            # ECOLIN_05405 and ECOLIN_01620 in AZ20
            SimpleLocation(ExactPosition(3888888), ExactPosition(3889127), strand=-1),
            SimpleLocation(ExactPosition(3888919), ExactPosition(3889127), strand=-1),
        ),
        'same_stop': (
            features['1-0006']['Rv2879c']['ratt'].location,
            features['1-0006']['Rv2880c']['ratt'].location,
        ),
        'inframe_pseudo_same_start_different_stop': (
            # 1-0006 RATT \pseudo SpmT with an internal stop codon
            SimpleLocation(ExactPosition(989728), ExactPosition(991202), strand=1),
            SimpleLocation(ExactPosition(989728), ExactPosition(990112), strand=1),
        ),
        'inframe_pseudo_different_start_different_stop': (
            # 1-0006 Rv2084 different start positions, but in-frame and there's an internal stop
            SimpleLocation(ExactPosition(2345919), ExactPosition(2346942), strand=1),
            SimpleLocation(ExactPosition(2345913), ExactPosition(2346915), strand=1),
        ),
        # 1-0006 combined CDS L_01053 and L_01054 compared to start-corrected
        'inframe_pseudo_on_nonoverlapping_part': (
            SimpleLocation(ExactPosition(1104310), ExactPosition(1105056), strand=1),
            SimpleLocation(ExactPosition(1105048), ExactPosition(1107616), strand=1),
        ),
        'overlapping_out_of_frame': (
            # 1-0006 espH and eccA1
            SimpleLocation(ExactPosition(4350755), ExactPosition(4351307), strand=1),
            # not a multiple of 3
            SimpleLocation(ExactPosition(4351299), ExactPosition(4353020), strand=1),
        ),
        'different_strand': (
            SimpleLocation(ExactPosition(0), ExactPosition(300), strand=1),
            SimpleLocation(ExactPosition(0), ExactPosition(300), strand=-1),
        ),
        'non_overlapping': (
            SimpleLocation(ExactPosition(0), ExactPosition(300), strand=1),
            SimpleLocation(ExactPosition(304), ExactPosition(335), strand=-1),
        ),
        'one_bp_apart': (
            # 1-0006's Rv0138 and Rv0139
            SimpleLocation(ExactPosition(163764), ExactPosition(164268), strand=1),
            SimpleLocation(ExactPosition(164268), ExactPosition(165291), strand=1),
        ),
        'compound_interval': (
            #PGAP annotation of 1-0013
            CompoundLocation([SimpleLocation(ExactPosition(1931), ExactPosition(2220), strand=1),
                              SimpleLocation(ExactPosition(2219), ExactPosition(3193), strand=1)], 'join'),
            CompoundLocation([SimpleLocation(ExactPosition(1931), ExactPosition(2258), strand=1),
                              SimpleLocation(ExactPosition(2207), ExactPosition(3193), strand=1)], 'join'),
        ),
    }
    expected = {
        'two_pseudos_same_start': True,
        'same_stop': True,
        'inframe_pseudo_same_start_different_stop': True,
        'inframe_pseudo_different_start_different_stop': True,
        'inframe_pseudo_on_nonoverlapping_part': False,
        'overlapping_out_of_frame': False,
        'different_strand': False,
        'non_overlapping': False,
        'one_bp_apart': False,
        'compound_interval': True,
    }

    assert compare.overlap_inframe(pairs[pair][0], pairs[pair][1]) == expected[pair]

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
            # nonconfl in 1
            [],
            # nonconfl in 2
            [
                inputs[case][1][1],
            ],
            # unique in 1
            [],
            # unique in 2
            [],
        ),
    }

    # comparison results should be identical regardless of which is called genome 1
    if swapped:
        expected[case] = list(expected[case])
        expected[case][0] = [(y, x) for x,y in expected[case][0]]
        expected[case][1] = [(y, x) for x,y in expected[case][1]]
        expected[case][2], expected[case][3] = expected[case][3], expected[case][2]
        expected[case][4], expected[case][5] = expected[case][5], expected[case][4]
        expected[case] = tuple(expected[case])

    assert compare.compare(f1_list, f2_list)[0:6] == expected[case]
