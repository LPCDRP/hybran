from collections import defaultdict, OrderedDict
from copy import deepcopy
import os

from Bio.SeqFeature import SeqFeature, SimpleLocation, FeatureLocation, ExactPosition, CompoundLocation

import pytest

from hybran import annomerge
from hybran import config
from hybran.bio import SeqIO
from hybran.util import keydefaultdict

from .data_features import *



def test_ref_fuse():
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features['H37Rv'])
    ref_id = 'H37Rv.NC_000962.3'

    assert annomerge.ref_annotation[f'{ref_id}::{ref_id}@@@PE_PGRS50::PE_PGRS49'].location == CompoundLocation([
        SimpleLocation(ExactPosition(3738157), ExactPosition(3742774), strand=-1, ref=ref_id),
        SimpleLocation(ExactPosition(3736983), ExactPosition(3738000), strand=-1, ref=ref_id)
    ], 'join')


@pytest.mark.parametrize('location', [
    'internal',
    'internal_minus',
    'beginning',
    'end_minus',
# not implemented
#    'beginning_linear',
#    'end_minus_linear',
])
@pytest.mark.skip("unused")
def test_upstream_context(location):
    # sequence generated with `makenucseq`
    seq = 'cgctaccgggtgcttgacta'
    n = 7
    locations = {
        'internal': FeatureLocation(ExactPosition(15),ExactPosition(17), strand=1),
        'internal_minus': FeatureLocation(ExactPosition(8),ExactPosition(10), strand=-1),
        'beginning': FeatureLocation(ExactPosition(5),ExactPosition(7), strand=1),
        'beginning_linear': FeatureLocation(ExactPosition(5),ExactPosition(7), strand=1),
        'end_minus': FeatureLocation(ExactPosition(15),ExactPosition(17), strand=-1),
        'end_minus_linear': FeatureLocation(ExactPosition(15),ExactPosition(17), strand=-1),
    }
    if location.endswith("_linear"):
        circular = False
    else:
        circular = True

    if circular:
        expected = {
            'internal': 'ggtgctt',
            'internal_minus': 'tgcttga',
            'beginning': 'tacgcta',
            'end_minus': 'ctacgct',
        }
    else:
        expected = {
            'beginning_linear': 'cgcta',
            'end_minus_linear': 'cta',
        }

    assert annomerge.upstream_context(
        locations[location],
        source_seq=seq,
        n=n,
        circular=circular
    ) == expected[location]

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
            FeatureLocation(ExactPosition(3888888), ExactPosition(3889127), strand=-1),
            FeatureLocation(ExactPosition(3888919), ExactPosition(3889127), strand=-1),
        ),
        'same_stop': (
            features['1-0006']['Rv2879c']['ratt'].location,
            features['1-0006']['Rv2880c']['ratt'].location,
        ),
        'inframe_pseudo_same_start_different_stop': (
            # 1-0006 RATT \pseudo SpmT with an internal stop codon
            FeatureLocation(ExactPosition(989728), ExactPosition(991202), strand=1),
            FeatureLocation(ExactPosition(989728), ExactPosition(990112), strand=1),
        ),
        'inframe_pseudo_different_start_different_stop': (
            # 1-0006 Rv2084 different start positions, but in-frame and there's an internal stop
            FeatureLocation(ExactPosition(2345919), ExactPosition(2346942), strand=1),
            FeatureLocation(ExactPosition(2345913), ExactPosition(2346915), strand=1),
        ),
        # 1-0006 combined CDS L_01053 and L_01054 compared to start-corrected
        'inframe_pseudo_on_nonoverlapping_part': (
            FeatureLocation(ExactPosition(1104310), ExactPosition(1105056), strand=1),
            FeatureLocation(ExactPosition(1105048), ExactPosition(1107616), strand=1),
        ),
        'overlapping_out_of_frame': (
            # 1-0006 espH and eccA1
            FeatureLocation(ExactPosition(4350755), ExactPosition(4351307), strand=1),
            # not a multiple of 3
            FeatureLocation(ExactPosition(4351299), ExactPosition(4353020), strand=1),
        ),
        'different_strand': (
            FeatureLocation(ExactPosition(0), ExactPosition(300), strand=1),
            FeatureLocation(ExactPosition(0), ExactPosition(300), strand=-1),
        ),
        'non_overlapping': (
            FeatureLocation(ExactPosition(0), ExactPosition(300), strand=1),
            FeatureLocation(ExactPosition(304), ExactPosition(335), strand=-1),
        ),
        'one_bp_apart': (
            # 1-0006's Rv0138 and Rv0139
            FeatureLocation(ExactPosition(163764), ExactPosition(164268), strand=1),
            FeatureLocation(ExactPosition(164268), ExactPosition(165291), strand=1),
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

    assert annomerge.overlap_inframe(pairs[pair][0], pairs[pair][1]) == expected[pair]

@pytest.mark.parametrize('gene_list', [
    'complementary_fragments',
    'complementary_fragments_one_unnamed',
    'fails_final_coord_check_inframe_overlap',
#    'seemingly_complete_fragment',
    'independent_copies',
])
@pytest.mark.skipif(not os.path.isfile("data/1-0006.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/2-0031.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.fasta"), reason="test reference genome sequence not available")
def test_fissionfuser(gene_list, tmp_path):
    inputs = {
        # post coordinate correction, so they both have the same start position
        'complementary_fragments' : [
            SeqFeature(
                FeatureLocation(ExactPosition(2275540), ExactPosition(2277261), strand=-1), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02173'],
                    'pseudo': [''],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(2276448), ExactPosition(2277261), strand=-1), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02174'],
                    'pseudo': [''],
                }
            ),
        ],
        'complementary_fragments_one_unnamed': [
            SeqFeature(
                FeatureLocation(ExactPosition(1104310), ExactPosition(1104607), strand=1), type='CDS', qualifiers={
                    'locus_tag': ['L_01053'],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(1104558), ExactPosition(1105056), strand=1), type='CDS', qualifiers={
                    'gene': ['Rv0986'],
                    'locus_tag': ['L_01054'],
                    'pseudo': [''],
                }
            ),
        ],
        'fails_final_coord_check_inframe_overlap' : [
            SeqFeature(
                # pre-correction: FeatureLocation(ExactPosition(3522323), ExactPosition(3523235), strand=-1)
                FeatureLocation(ExactPosition(3522323), ExactPosition(3523418), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03351'],
                    'pseudo': [''],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(3523256), ExactPosition(3523418), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03352'],
                    'pseudo': [''],
                }
            ),
        ],
        'seemingly_complete_fragment': [
        ],
        'independent_copies': [
            SeqFeature(
                FeatureLocation(ExactPosition(1483334), ExactPosition(1484933), strand=-1), type='CDS', qualifiers={
                    'locus_tag': ['L_01415'],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(1484939), ExactPosition(1486616), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv1319c'],
                    'locus_tag': ['L_01416'],
                }
            ),
        ],
    }
    source_genomes = {
        'complementary_fragments':'1-0006',
        'complementary_fragments_one_unnamed':'1-0006',
        'fails_final_coord_check_inframe_overlap':'2-0031',
        'independent_copies':'1-0006',
    }
    config.hybran_tmp_dir = tmp_path
    ref_genome = defaultdict(lambda :'H37Rv')
    source_genome = source_genomes[gene_list]
    record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
    annomerge.genetic_code = 11
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[gene_list]])

    for f in inputs[gene_list]:
        f.source = 'H37Rv.NC_000962.3' # TODO - clean this up
        f.ref = record_sequence.id
        f.references = {record_sequence.id: record_sequence.seq}

    expected = {
        'complementary_fragments': ([
            SeqFeature(
                FeatureLocation(ExactPosition(2275540), ExactPosition(2277261), strand=-1, ref='1'), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02174'],
                    'note': [
                        ('Hybran/Pseudoscan: Internal stop detected in the following '
                         'codon(s): 270 283 353 380 | Locus has invalid reading frame-- not '
                         'divisible by three | Locus has reference-corresponding start and '
                         'end | Poor blastp match at 95% identity and 95% coverage thresholds '
                         '| Locus is 1 base pair(s) shorter than the reference')
                    ],
                    'pseudo': [''],
                }
            )],
            [{
                'feature':inputs['complementary_fragments'][0],
                'evid':'overlapping_inframe',
                'remark':'combined fission fragments',
            }],
        ),
        'complementary_fragments_one_unnamed': ([
            SeqFeature(
                FeatureLocation(ExactPosition(1104310), ExactPosition(1105056), strand=1, ref='1'), type='CDS', qualifiers={
                    'gene': ['Rv0986'],
                    'locus_tag': ['L_01054'],
                    'note': [
                        ('Hybran/Pseudoscan: Internal stop detected in the following '
                         'codon(s): 98 168 199 208 215 223 227 | Locus has invalid reading '
                         'frame-- not divisible by three | Locus has reference-corresponding '
                         'start and end | Poor blastp match at 95% identity and 95% coverage '
                         'thresholds | Locus is 1 base pair(s) shorter than the reference')
                    ],
                    'pseudo': [''],
                }
            )],
            [{
                'feature':inputs['complementary_fragments_one_unnamed'][0],
                'evid':'complementary_fragments',
                'remark':'combined fission fragments',
            }],
        ),
        'fails_final_coord_check_inframe_overlap': ([
            SeqFeature(
                FeatureLocation(ExactPosition(3522323), ExactPosition(3523418), strand=-1, ref='1'), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03352'],
                    'note': [
                        ('Hybran/Pseudoscan: Internal stop detected in the following '
                         'codon(s): 53 | Locus has valid reading frame | Locus does not have '
                         'reference-corresponding end | Poor blastp match at 95% identity and '
                         '95% coverage thresholds | Locus has a delayed stop codon')
                    ],
                    'pseudo': [''],
                }
            )],
            [{
                'feature':inputs['fails_final_coord_check_inframe_overlap'][0],
                'evid':'overlapping_inframe',
                'remark':"combined fission fragments",
            }],
        ),
        'independent_copies': (
            inputs['independent_copies'],
            [],
        ),
    }

    # set the rival feature as the passing one. for these test cases, we're always looking at pairs
    if expected[gene_list][1]:
        expected[gene_list][1][0]['superior'] = expected[gene_list][0][0]

    assert annomerge.fissionfuser(
        inputs[gene_list],
        seq_ident=95,
        seq_covg=95,
    ) == expected[gene_list]

@pytest.mark.parametrize('gene_list', [
    'misannotation_false_delayed_stop',
    'redundant_double_hybrid_fusion',
    'misannotation_both_nonpseudo',
    'misannotation_one_pseudo',
    'idempotence',
])
def test_fusionfisher(gene_list):
    source_genomes = {
        'misannotation_false_delayed_stop':'1-0006',
        'redundant_double_hybrid_fusion':'AZ20',
        'misannotation_both_nonpseudo': 'AZ20',
        'misannotation_one_pseudo': 'AZ20',
        'idempotence': 'AZ20',
    }
    # create a dummy feature dictionary for the test case that is not in focus in the current invocation to avoid a KeyError when loading all expected inputs and results
    source_features = defaultdict(lambda :defaultdict(dict))
    source_features.update(features[source_genomes[gene_list]])
    inputs = {
        'misannotation_false_delayed_stop': [
            source_features['Rv0074']['ratt'],
            source_features['Rv0071']['ratt'],
        ],
        'redundant_double_hybrid_fusion': [
            source_features['AZ20_03933']['ratt'],
            source_features['AZ20_03933']['prokka'],
        ],
        'misannotation_both_nonpseudo': [
            source_features['ECOLIN_18975']['ratt'],
            source_features['ECOLIN_18965']['ratt'],
        ],
        'misannotation_one_pseudo': [
            source_features['ECOLIN_25305']['ratt'],
            source_features['garD']['ratt'],
        ],
        'idempotence': [
            deepcopy(source_features['ECOLIN_24700']['ratt']),
            deepcopy(source_features['ECOLIN_12095']['prokka']),
        ],
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome.update({
        'redundant_double_hybrid_fusion': 'nissle-hybrid',
        'misannotation_both_nonpseudo': 'nissle-hybrid',
        'misannotation_one_pseudo': 'nissle-hybrid',
        'idempotence': 'nissle-hybrid',
    })
    source_genome = source_genomes[gene_list]

    record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
    for f in inputs[gene_list]:
        f.references = {record_sequence.id: record_sequence.seq}
        for part in f.location.parts:
            part.ref = record_sequence.id
    annomerge.genetic_code = 11
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[gene_list]])
    expected = {
        'misannotation_false_delayed_stop': (
            [ source_features['Rv0074']['ratt'] ],
            [ ],
            [{
                'feature':source_features['Rv0071']['ratt'],
                'evid':'putative_misannotation',
                'remark':"Has no reference-corresponding stop, while rival feature does, and both share the same stop position.",
            }],
        ),
        'redundant_double_hybrid_fusion': (
            [ source_features['AZ20_03933']['ratt'] ],
            [ ],
            [{
                'feature':source_features['AZ20_03933']['prokka'],
                'evid':'redundant_fusion_member',
                'remark':"Same locus as rival (fusion) gene and name already included as fusion element there.",
            }],
        ),
        'misannotation_both_nonpseudo': (
            [ source_features['ECOLIN_18975']['ratt'] ],
            [ ],
            [{
                'feature':source_features['ECOLIN_18965']['ratt'],
                'evid':'shorter',
                'remark':"Both have reference-corresponding start codons and have reference-corresponding stop codons.",
            }],
        ),
        'misannotation_one_pseudo': (
            [source_features['ECOLIN_25305']['ratt']],
            [ ],
            [{
                'feature':source_features['garD']['ratt'],
                'evid':'putative_misannotation',
                'remark':"Has no reference-corresponding coordinates, while rival feature has a reference-corresponding start, and both share the same stop position.",
            }],
        ),
        'idempotence': (
            deepcopy(inputs['idempotence']),
            [ inputs['idempotence'][0] ],
            [ ],
        ),
    }
    if gene_list == 'idempotence':
        expected['idempotence'][0][1].qualifiers['note'] = [
            "Upstream gene ECOLIN_24700|ECOLIN_24700::ECOLIN_12095 conjoins with this one."
        ]

    # set the rival feature as the passing one. for these test cases, we're always looking at pairs
    if expected[gene_list][2]:
        expected[gene_list][2][0]['superior'] = expected[gene_list][0][0]

    assert annomerge.fusionfisher(
        deepcopy(inputs[gene_list]),
    ) == expected[gene_list]


@pytest.mark.skip("superseded by fusionfisher")
def test_identify_conjoined_genes():
    ratt_features = [
        features['1-0006']['PPE5']['ratt'],
        features['1-0006']['PPE6']['ratt'],
        features['1-0006']['Rv2879c']['ratt'],
        features['1-0006']['Rv2880c']['ratt'],
    ]

    # compare locations because SeqFeature objects can't be directly compared properly.
    # https://github.com/biopython/biopython/issues/3874
    expected_locs = [
        features['1-0006']['PPE6']['ratt'].location
    ]

    observed = annomerge.identify_conjoined_genes(ratt_features)

    assert [_.location for _ in observed] == expected_locs


@pytest.mark.skip("superseded by fusionfisher")
def test_identify_merged_genes():
    ratt_features = [
        features['1-0006']['PPE5']['ratt'],
        features['1-0006']['PPE6']['ratt'],
        features['1-0006']['Rv2879c']['ratt'],
        features['1-0006']['Rv2880c']['ratt'],
    ]

    expected = {
        -1:[
            (int(ratt_features[2].location.start) , int(ratt_features[2].location.end)),
        ],
        1:[
        ]
    }

    assert annomerge.identify_merged_genes(ratt_features) == (expected, True)

def test_liftover_annotation():
    abinit = SeqFeature(
        FeatureLocation(ExactPosition(0), ExactPosition(10), strand=1), type='CDS',
        qualifiers={
            'locus_tag' : ['L_00017'],
            'gene' : ['L_00017'],
            'inference' : [
                'ab initio prediction:Prodigal:002006',
                'similar to AA sequence:UniProtKB:Q9WYC4'
            ],
            'codon_start' : ['1'],
            'transl_table': ['4'],
            'product': ['putative ABC transporter ATP-binding protein'],
            'protein_id': ['C:L_00017'],
            'db_xref': ['COG:COG1132'],
            'translation': ['MQL']
        }
    )
    ref = SeqFeature(
        FeatureLocation(ExactPosition(0), ExactPosition(10), strand=1), type='CDS',
        qualifiers={
            'locus_tag': ['MG_RS00085'],
            'gene': ['MG_RS00085'],
            'old_locus_tag': ['MG_015'],
            'note': [
                '*inference: COORDINATES: similar to AA sequence:RefSeq:WP_014894097.1',
                'Derived by automated computational analysis using gene prediction method: Protein Homology.',
            ],
            'codon_start': ['1'],
            'transl_table': ['4'],
            'product': ['ABC transporter ATP-binding protein/permease'],
            'protein_id': ['WP_010869291.1'],
            'translation': ['MEG'],
        }
    )
    expected = {
        'locus_tag' : ['L_00017'],
        'gene' : ['MG_RS00085'],
        'inference' : [
            'ab initio prediction:Prodigal:002006',
            'alignment:blastp'
        ],
        'note': [
            '*inference: COORDINATES: similar to AA sequence:RefSeq:WP_014894097.1',
            'Derived by automated computational analysis using gene prediction method: Protein Homology.',
        ],
        'codon_start' : ['1'],
        'transl_table': ['4'],
        'product': ['ABC transporter ATP-binding protein/permease'],
        'protein_id': ['WP_010869291.1'],
        'db_xref': ['COG:COG1132'],
        'translation': ['MQL']
    }
    annomerge.liftover_annotation(abinit, ref, inference='alignment:blastp')
    assert abinit.qualifiers == expected

@pytest.mark.parametrize('feature_type,fix_start,fix_stop,seek_stop', [
    ['abinit_start_bad_minus', True, False, None],
    ['bad_translation', True, False, None],
    ['tricky_found_low', True, False, None],
    ['good_start_stop_deletion', True, False, None],
    ['bad_mismatch_check', True, True, None],
    ['bad_mismatch_check2', True, True, None],
    ['ref_start_frameshift', True, False, None],
    ['bad_start_good_padding', True, False, None],
    ['ratt_pseudo_pgrs', True, False, None],
    ['same_start_alt_stop_1', False, False, None],
    ['same_start_alt_stop_2', False, False, None],
    ['same_start_alt_stop_2_fix', True, True, False],
    ['bad_start_stop_nofix_pseudo', False, False, None],
    ['bad_start_stop_fix_pseudo', True, False, None],
    ['inverted_join_ecoli', True, True, None],
    ['gene_fusion', True, True, None],
    ['end_greater_than_start', True, True, None],
    ['extend_gap_penalty_delay_stop', True, True, None],
    ['insertion_final_interval', True, True, None],
    ['early_del_altered_scoring', False, False, None],
    ['alt_start_delicate_scoring', True, False, None],
    ['first_codon_snp_valid_start', True, True, None],
])
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.gbk"), reason="test reference annotation not available")
@pytest.mark.skipif(not os.path.isfile("data/nissle-hybrid.gbk"), reason="test reference annotation not available")
@pytest.mark.skipif(not os.path.isfile("data/PAO1_107.gbk"), reason="test reference annotation not available")
def test_coord_check(feature_type, fix_start, fix_stop, seek_stop):
    #prokka for Rv2300c and Rv3181c
    source_genome = {
        'abinit_start_bad_minus':'1-0006',
        'bad_translation':'1-0006',
        'tricky_found_low':'1-0006',
        'good_start_stop_deletion':'1-0006',
        'bad_mismatch_check':'1-0006',
        'bad_mismatch_check2':'1-0006',
        'ref_start_frameshift':'1-0006',
        'bad_start_good_padding':'2-0031',
        'ratt_pseudo_pgrs':'1-0006',
        'same_start_alt_stop_1':'1-0006',
        'same_start_alt_stop_2':'1-0006',
        'same_start_alt_stop_2_fix':'1-0006',
        'bad_start_stop_nofix_pseudo':'1-0006',
        'bad_start_stop_fix_pseudo':'1-0006',
        'inverted_join_ecoli':'AZ20',
        'gene_fusion':'1-0006',
        'end_greater_than_start':'SEA08151',
        'extend_gap_penalty_delay_stop':'1-0006',
        'insertion_final_interval':'PAK',
        'early_del_altered_scoring':'1-0006',
        'alt_start_delicate_scoring':'1-0006',
        'first_codon_snp_valid_start':'1-0006',

    }
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome.update({
        'inverted_join_ecoli': 'nissle-hybrid',
        'insertion_final_interval': 'PAO1_107',
    })

    test_features = {
        'abinit_start_bad_minus': features[source_genome['abinit_start_bad_minus']]['Rv3181c']['abinit'],
        'bad_translation': features[source_genome['bad_translation']]['Rv3777']['abinit'],
        'tricky_found_low': features[source_genome['tricky_found_low']]['PPE38']['abinit'],
        'good_start_stop_deletion': features[source_genome['good_start_stop_deletion']]['Rv3785']['abinit'],
        'bad_mismatch_check': features[source_genome['bad_mismatch_check']]['Rv1877']['ratt'],
        'bad_mismatch_check2': features[source_genome['bad_mismatch_check2']]['Rv3327']['ratt'],
        'ref_start_frameshift' : features[source_genome['ref_start_frameshift']]['PPE47']['abinit'],
        'bad_start_good_padding' : features[source_genome['bad_start_good_padding']]['PPE34']['abinit'],
        'ratt_pseudo_pgrs': features[source_genome['ratt_pseudo_pgrs']]['PE_PGRS50']['ratt_raw'],
        'same_start_alt_stop_1': features[source_genome['same_start_alt_stop_1']]['Rv2879c']['ratt'],
        'same_start_alt_stop_2': features[source_genome['same_start_alt_stop_2']]['Rv2880c']['ratt'],
        'same_start_alt_stop_2_fix': features[source_genome['same_start_alt_stop_2']]['Rv2880c']['ratt'],
        'bad_start_stop_nofix_pseudo': features[source_genome['bad_start_stop_nofix_pseudo']]['PE10']['ratt'],
        'bad_start_stop_fix_pseudo': features[source_genome['bad_start_stop_nofix_pseudo']]['PE10']['ratt'],
        'inverted_join_ecoli': features[source_genome['inverted_join_ecoli']]['secD']['ratt'],
        'gene_fusion': features[source_genome['gene_fusion']]['PE_PGRS50']['final'],
        'end_greater_than_start': features[source_genome['end_greater_than_start']]['lpqG']['ratt'],
        'extend_gap_penalty_delay_stop': features[source_genome['extend_gap_penalty_delay_stop']]['Rv0325']['ratt'],
        'insertion_final_interval': features[source_genome['insertion_final_interval']]['PA2452']['ratt'],
        'early_del_altered_scoring': features[source_genome['early_del_altered_scoring']]['accE5']['ratt'],
        'alt_start_delicate_scoring': features[source_genome['alt_start_delicate_scoring']]['Rv3611']['ratt'],
        'first_codon_snp_valid_start': features[source_genome['first_codon_snp_valid_start']]['Rv2023A']['ratt'],

    }

    record_sequence = list(SeqIO.parse(f'data/{source_genome[feature_type]}.fasta', 'fasta'))[0]
    test_features[feature_type].ref = record_sequence.id
    test_features[feature_type].references = {record_sequence.id: record_sequence.seq}
    annomerge.genetic_code = 11
    annomerge.corrected_orf_report = []
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[feature_type]])

    feature = test_features[feature_type]
    ref_feature = annomerge.ref_annotation[
        annomerge.key_ref_gene(test_features[feature_type].source, test_features[feature_type].qualifiers['gene'][0])
    ]

    expected = {
        'abinit_start_bad_minus': {
            'results':[(True, True), FeatureLocation(3548089, 3548542, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':False,
        },
        'bad_translation': {
            'results':[(True, True), FeatureLocation(4230767, 4231754, strand=1, ref='1')],
            'og_de':False,
            'corr_de':False,
        },
        'tricky_found_low': {
            'results':[(False, True), FeatureLocation(2636045, 2637140, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'good_start_stop_deletion': {
            'results':[(True, True), FeatureLocation(4239393, 4240378, strand=1, ref='1')],
            'og_de':False,
            'corr_de':False,
        },
        'bad_mismatch_check': {
            'results':[(True, False), FeatureLocation(2112334, 2114834, strand=1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'bad_mismatch_check2': {
            'results':[(False, False), FeatureLocation(3707086, 3709176, strand=1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'ref_start_frameshift': {
            'results':[(True, True), FeatureLocation(3374234, 3375312, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':False,
        },
        'bad_start_good_padding': {
            'results':[(False, True), FeatureLocation(1927378, 1927894, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'ratt_pseudo_pgrs': {
            'results':[(True,False), FeatureLocation(3741108, 3746955, strand=-1, ref='1')],
            'og_de':True,
            'corr_de':None,
        },
        'same_start_alt_stop_1': {
            'results':[(False, True), FeatureLocation(3182302, 3183397, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'same_start_alt_stop_2': {
            'results':[(True, False), FeatureLocation(3182302, 3183397, strand=-1, ref='1')],
            'og_de':True,
            'corr_de':None,
        },
        'same_start_alt_stop_2_fix': {
            'results':[(True, True), FeatureLocation(3182570, 3183397, strand=-1, ref='1')],
            'og_de':True,
            'corr_de':False,
        },
        'bad_start_stop_nofix_pseudo': {
            'results':[(False, False), FeatureLocation(1217413, 1217872, strand=1, ref='1')],
            'og_de':True,
            'corr_de':None,
        },
        'bad_start_stop_fix_pseudo':  {
            'results':[(True, False), FeatureLocation(1217428, 1217872, strand=1, ref='1')],
            'og_de':True,
            'corr_de':True,
        },
        'inverted_join_ecoli':  {
            'results':[(False, False), FeatureLocation(3714209, 3716770, strand=-1, ref='1')],
            'og_de':True,
            'corr_de':None,
        },
        'gene_fusion':  {
            'results':[(True, True), FeatureLocation(3741108, 3746955, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'end_greater_than_start':  {
            'results':[(False, False), FeatureLocation(4064133, 4064322, strand=1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'extend_gap_penalty_delay_stop':  {
            'results':[(True, False), FeatureLocation(392626, 393316, strand=1, ref='1')],
            'og_de':True,
            'corr_de':None,
        },
        'insertion_final_interval':  {
            'results':[(True, False), FeatureLocation(2799744, 2801325, strand=1, ref='refseq|NZ_LR657304.1|chromosome1')],
            'og_de':True,
            'corr_de':None,
        },
        'early_del_altered_scoring':  {
            'results':[(True, True), FeatureLocation(3659371, 3659770, strand=1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'alt_start_delicate_scoring':  {
            'results':[(True, True), FeatureLocation(4060333, 4061209, strand=1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
        'first_codon_snp_valid_start': {
            'results':[(True, True), FeatureLocation(2264663, 2265122, strand=-1, ref='1')],
            'og_de':False,
            'corr_de':None,
        },
    }
    results = annomerge.coord_check(feature, ref_feature, fix_start=fix_start, fix_stop=fix_stop, seek_stop=seek_stop)
    assert [
        [results, feature.location],
        feature.og.de, feature.corr.de
    ] == [
        expected[feature_type]['results'],
        expected[feature_type]['og_de'],
        expected[feature_type]['corr_de'],
    ]

@pytest.mark.parametrize('feature_type, seq_ident, seq_covg, attempt_rescue', [
    ['small_badstop_fix_pseudo', 95, 95, True],
    ['small_badstart_fix_nopseudo', 95, 95, True],
    ['small_badstop_fix_pseudo_frameshift', 95, 95, True],
    ['frameshift_alt_stop_pseudo', 95, 95, True],
    ['good_start_stop_frameshift_pseudo', 95, 95, True],
    ['bad_start_stop_nofix_pseudo', 95, 95, True],
    ['fix_stop_valid_broken_stop', 95, 95, True],
    ['sensitive_padding_fix_pseudo', 95, 95, True],
    ['deletion_in_middle_fix_pseudo', 95, 95, True],
    ['fix_start_stop_nonpseudo', 95, 95, True],
    ['good_start_stop_fix_pseudo', 95, 95, True],
    ['inframe_deletion_in_middle', 95, 95, True],
    ['good_blast_still_repairable', 95, 95, True],
    ['start_correction_induces_delayed_stop', 95, 95, True],
    ['start_correction_induces_delayed_stop2', 95, 95, True],
    ['reject_coord_correction', 95, 95, True],
])
def test_pseudoscan(feature_type, seq_ident, seq_covg, attempt_rescue, tmp_path):
    ref_genome = defaultdict(lambda :'H37Rv')
    source_genome = {
        'small_badstop_fix_pseudo':'1-0006',
        'small_badstart_fix_nopseudo':'1-0006',
        'small_badstop_fix_pseudo_frameshift':'1-0006',
        'frameshift_alt_stop_pseudo': '1-0006',
        'good_start_stop_frameshift_pseudo': '1-0006',
        'bad_start_stop_nofix_pseudo': '1-0006',
        'fix_stop_valid_broken_stop': '1-0006',
        'sensitive_padding_fix_pseudo': '1-0006',
        'deletion_in_middle_fix_pseudo': '1-0006',
        'fix_start_stop_nonpseudo': '1-0006',
        'good_start_stop_fix_pseudo': '1-0006',
        'inframe_deletion_in_middle': '1-0009',
        'good_blast_still_repairable': '1-0006',
        'start_correction_induces_delayed_stop': '1-0006',
        'start_correction_induces_delayed_stop2': '1-0006',
        'reject_coord_correction': '1-0006',
    }

    test_features = {
        'small_badstop_fix_pseudo': features[source_genome['small_badstop_fix_pseudo']]['Rv0061c']['ratt'],
        'small_badstart_fix_nopseudo': features[source_genome['small_badstart_fix_nopseudo']]['galTb']['ratt'],
        'small_badstop_fix_pseudo_frameshift':features[source_genome['small_badstop_fix_pseudo_frameshift']]['Rv1075c']['ratt'],
        'frameshift_alt_stop_pseudo': features[source_genome['frameshift_alt_stop_pseudo']]['mce1R']['ratt'],
        'good_start_stop_frameshift_pseudo': features[source_genome['good_start_stop_frameshift_pseudo']]['PPE6']['ratt'],
        'bad_start_stop_nofix_pseudo': features[source_genome['bad_start_stop_nofix_pseudo']]['PE10']['ratt'],
        'fix_stop_valid_broken_stop': features[source_genome['fix_stop_valid_broken_stop']]['Rv2079']['ratt'],
        'sensitive_padding_fix_pseudo': features[source_genome['sensitive_padding_fix_pseudo']]['Rv2081c']['ratt'],
        'deletion_in_middle_fix_pseudo': features[source_genome['deletion_in_middle_fix_pseudo']]['Rv2134c']['ratt'],
        'fix_start_stop_nonpseudo': features[source_genome['fix_start_stop_nonpseudo']]['Rv2280']['ratt'],
        'good_start_stop_fix_pseudo': features[source_genome['good_start_stop_fix_pseudo']]['Rv2437']['ratt'],
        'inframe_deletion_in_middle': features[source_genome['inframe_deletion_in_middle']]['PPE54']['ratt'],
        'good_blast_still_repairable': features[source_genome['good_blast_still_repairable']]['dnaA']['abinit'],
        'start_correction_induces_delayed_stop': features[source_genome['start_correction_induces_delayed_stop']]['Rv1225c']['abinit_raw'],
        'start_correction_induces_delayed_stop2': features[source_genome['start_correction_induces_delayed_stop2']]['Rv2561']['ratt_raw'],
        'reject_coord_correction': features[source_genome['reject_coord_correction']]['pks6']['abinit'],
    }

    feature = test_features[feature_type]
    config.hybran_tmp_dir = tmp_path
    ref_feature = ref_features[ref_genome[feature_type]][
        annomerge.key_ref_gene(test_features[feature_type].source, test_features[feature_type].qualifiers['gene'][0])
    ]
    record_sequence = list(SeqIO.parse(f'data/{source_genome[feature_type]}.fasta', 'fasta'))[0]
    test_features[feature_type].ref = record_sequence.id
    test_features[feature_type].references = {record_sequence.id: record_sequence.seq}
    annomerge.genetic_code = 11
    annomerge.corrected_orf_report = []
    expected = {
        'small_badstop_fix_pseudo': [True, FeatureLocation(66693, 67032, strand=-1, ref='1')],
        'small_badstart_fix_nopseudo': [False, FeatureLocation(712762, 713308, strand=1, ref='1')],
        'small_badstop_fix_pseudo_frameshift': [True, FeatureLocation(1202098, 1203042, strand=-1, ref='1')],
        'frameshift_alt_stop_pseudo' : [True, FeatureLocation(192566, 193259, strand=-1, ref='1')],
        'good_start_stop_frameshift_pseudo' : [True, FeatureLocation(366753, 376299, strand=-1, ref='1')],
        'bad_start_stop_nofix_pseudo': [True, FeatureLocation(1217428, 1217872, strand=1, ref='1')],
        'fix_stop_valid_broken_stop': [True, FeatureLocation(2339465, 2341436, strand=1, ref='1')],
        'sensitive_padding_fix_pseudo': [True, FeatureLocation(2342175, 2342617, strand=-1, ref='1')],
        'deletion_in_middle_fix_pseudo': [True, FeatureLocation(2397532, 2398106, strand=-1, ref='1')],
        'fix_start_stop_nonpseudo': [False, FeatureLocation(2552370, 2553750, strand=1, ref='1')],
        'good_start_stop_fix_pseudo': [True, FeatureLocation(2735891, 2736310, strand=1, ref='1')],
        'inframe_deletion_in_middle': [False, FeatureLocation(3730264, 3741334, strand=-1, ref='1')],
        'good_blast_still_repairable': [False, FeatureLocation(0, 1524, strand=1, ref='1')],
        'start_correction_induces_delayed_stop': [True, FeatureLocation(1370377, 1371385, strand=-1, ref='1')],
        'start_correction_induces_delayed_stop2': [True, FeatureLocation(2881559, 2882297, strand=1, ref='1')],
        'reject_coord_correction': [False, FeatureLocation(486288, 490458, strand=1, ref='1')],
    }
    results = annomerge.pseudoscan(feature, ref_feature, seq_ident, seq_covg, attempt_rescue)
    assert [results, feature.location] == expected[feature_type]

def test_populate_gaps():
    intergenic_positions = [
        (1525, 3409, '+'),
        #(4619, 4637, '+'),
        (6356, 6597, '+'),
        # (8626, 8659, '+'),
        # (11177, 11271, '+'),
        (12187, 13843, '+'),
        # (13670, 14490, '-'),
        (14375, 15446, '+'),
        (15354, 16947, '-'),
        (16236, 16271, '+'),
        (160607, 160608, '+'),
        (400113, 400779, '+'),
        # (4416974, 4417044, '-'),
    ]

    ratt_pre_intergene = {
        (1524, 1): SeqFeature(
            FeatureLocation(0, 1524, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0001'],'gene':['dnaA']}
        ),
        (6355, 1): SeqFeature(
            FeatureLocation(5791, 6355, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0004'],'gene':['Rv0004']}
        ),
        (12186, 1): SeqFeature(
            FeatureLocation(11271, 12186, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0007'],'gene':['Rv0007']}
        ),
        (16235, 1): SeqFeature(
            FeatureLocation(15447, 16235, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0012'],'gene':['Rv0012']}
        ),
        (160606, 1): SeqFeature(
            FeatureLocation(160213, 160606, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0134'],'gene':['ephF'],'pseudo':['']}
        ),
        (400112, 1): SeqFeature(
            FeatureLocation(399245, 400112, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0334'],'gene':['rmlA']}
        )
    }
    ratt_post_intergene = {
        (3409, 1): SeqFeature(FeatureLocation(3409, 4618, strand=1), type='CDS'),
        (6597, 1): SeqFeature(
            FeatureLocation(6597, 8625, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0005'],'gene':['gyrB']}
        ),
        (16272, 1): SeqFeature(
            FeatureLocation(16272, 16970, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0013'],'gene':['trpG']}
        ),
        (160608, 1): SeqFeature(
            FeatureLocation(160608, 161115, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0134'],'gene':['ephF'],'pseudo':['']}
        ),
    }
    abinit_features = {
        (0, 4419608, 1): SeqFeature(FeatureLocation(0, 4419608, strand=1), type='source'),
        # overlaps before unannotated region
        (33, 1524, 1): SeqFeature(FeatureLocation(33, 1524, strand=1), type='CDS'),
        # squarely within unannotated region
        (1645, 1972, 1): SeqFeature(FeatureLocation(1645, 1972, strand=1), type='CDS'),
        # begins in unannotated region and ends past it
        (6480, 8625, 1): SeqFeature(FeatureLocation(6480, 8625, strand=1), type='CDS'),
        # overlaps before unannotated region, but not immediately before it (overlapping gene with different reading frame)
        (11320, 11518, 1): SeqFeature(FeatureLocation(11320, 11518, strand=1), type='CDS'),
        # ends at same position as a RATT gene, but starts in an intergenic region.
        # this test case is actually special because its intergenic region follows a minus strand region,
        (15491, 16235, 1): SeqFeature(FeatureLocation(15491, 16235, strand=1), type='CDS'),
        # bridges an unannotated region
        (160602, 161115, 1): SeqFeature(FeatureLocation(160602, 161115, strand=1), type='CDS'),
        # overlaps before unannotated region and into it
        (400074, 400686, 1): SeqFeature(FeatureLocation(400074, 400686, strand=1), type='CDS'),
    }

    expected_keepers = [
        abinit_features[(1645, 1972, 1)],
    ]
    expected_conflicts = defaultdict(list)
    # The [] + [] is to separate out the result that's expected only because I didn't specify every intermediate interval.
    # The second list would not be expected in the output if we were using a complete input dataset, but is correct
    # for what's provided here.
    expected_conflicts[(33, 1524, 1)] = [(0, 1524, 1)]
    expected_conflicts[(6480, 8625, 1)] = [(6597, 8625, 1)] + [(11271, 12186, 1)]
    expected_conflicts[(11320, 11518, 1)] = [(11271, 12186, 1)]
    expected_conflicts[(15491, 16235, 1)] = [(15447, 16235, 1)]
    expected_conflicts[(160602, 161115, 1)] = [(160213, 160606, 1), (160608, 161115, 1)] + [(399245, 400112, 1)]
    expected_conflicts[(400074, 400686, 1)] = [(399245, 400112, 1)]

    (keepers, conflicts) = annomerge.populate_gaps(
        abinit_features,
        intergenic_positions,
        ratt_pre_intergene,
        ratt_post_intergene,
    )
    # SeqFeatures can't be compared directly yet, so just compare the locations
    assert (([_.location for _ in keepers], conflicts)
            == ([_.location for _ in expected_keepers], expected_conflicts))

@pytest.mark.parametrize(
    "case",
    [
        'double_overlap',
        'overlapping_fused_ref',
    ]
)
def test_find_inframe_overlaps(case):
    # double-overlap
    # issues seen when annotating AZ20 using re-annotated NISSLE as reference
    nissleorf0212 = SeqFeature(
        FeatureLocation(2602635, 2603246, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['ECOLIN_26724'],'gene':['NISSLEORF0212']}
    )
    nissleorf0346 = SeqFeature(
        FeatureLocation(2603304, 2603475, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['ECOLIN_26723'],'gene':['NISSLEORF0346']}
    )
    abinit_nissleorf0212 = SeqFeature(
        FeatureLocation(2602635, 2603475, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['L_02518'],'gene':['NISSLEORF0212']}
    )

    # 1-0006. abinit overlapping fused reference.
    # not rejected here since names don't match. It's left to the inclusion criteria.
    abinit_ppe6 = SeqFeature(
        FeatureLocation(366753, 376299, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['L_00329']}
    )
    ppe5 = SeqFeature(
        FeatureLocation(366753, 373353, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['Rv0304c'], 'gene':['PPE5']}
    )
    ppe6 = SeqFeature(
        FeatureLocation(366753, 376299, strand=-1),
        type='CDS',
        qualifiers={'locus_tag':['Rv0305c'], 'gene':['PPE6']}
    )

    def loc_triplet(feature):
        return (feature.location.start, feature.location.end, feature.location.strand)
    def dictate(features):
        return {loc_triplet(_):_ for _ in features}

    ratt_features = {
        'double_overlap':[nissleorf0212, nissleorf0346],
        'overlapping_fused_ref': [ppe5, ppe6],
    }
    abinit_features = {
        'double_overlap':[abinit_nissleorf0212],
        'overlapping_fused_ref': [abinit_ppe6],
    }
    expected_keep = {
        'double_overlap': [abinit_nissleorf0212],
        'overlapping_fused_ref': [abinit_ppe6],
    }
    expected_rejects = {
        'double_overlap': [],
        'overlapping_fused_ref': [],
    }
    expected_overlaps = {
        'double_overlap': defaultdict(list, {loc_triplet(abinit_nissleorf0212): [loc_triplet(_) for _ in [nissleorf0212, nissleorf0346]]}),
        'overlapping_fused_ref': defaultdict(list, {loc_triplet(abinit_ppe6): [loc_triplet(_) for _ in [ppe5, ppe6]]}),
    }

    abinit_features_dictionary = dictate(abinit_features[case])

    assert annomerge.find_inframe_overlaps(ratt_features[case], abinit_features_dictionary) == \
        (dictate(expected_keep[case]), expected_overlaps[case], expected_rejects[case])

@pytest.mark.parametrize('pair', [
    'ratt_better',
    'ratt_better_coverage',
    'pseudo_vs_nonpseudo',
    'overlapping_unnamed',
    'abinit_better',
    # https://gitlab.com/LPCDRP/hybran/-/issues/57
    'overlapping_different_names_ratt_better',
    'overlapping_different_names_abinit_better',
    'ratt_join_vs_prokka_bad_start',
    'prokka_gene_fusion',
])
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.gbk"), reason="test reference annotation not available")
def test_check_inclusion_criteria(pair, tmp_path):
    source_genome = {
        'ratt_better':'1-0006',
        'ratt_better_coverage':'1-0006',
        'pseudo_vs_nonpseudo':'1-0006',
        'overlapping_unnamed':'1-0006',
        'abinit_better':'4-0041',
        'corresponding_non_cds':'4-0041',
        'overlapping_different_names_ratt_better':'1-0006',
        'overlapping_different_names_abinit_better':'1-0006',
        'ratt_join_vs_prokka_bad_start':'1-0006',
        'prokka_gene_fusion':'1-0006',
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    pairs = {
        'ratt_better': ('dnaA', 'dnaA'),
        'ratt_better_coverage': ('Rv1453', 'Rv1453'), # The RATT annotation's upstream context has no hit to the reference's despite being 100% identitical...
        'pseudo_vs_nonpseudo': ('Rv0007','Rv0007'),
        'overlapping_unnamed': ('dnaA', 'gyrB'),
        'corresponding_non_cds': ('rrf', 'rrf'),
        'abinit_better': ('Rv1718', 'Rv1718'),
        'overlapping_different_names_ratt_better': ('Rv1945', 'Rv1945'),
        'overlapping_different_names_abinit_better': ('Rv2180c', 'ORF0004'),
        'ratt_join_vs_prokka_bad_start':('pip', 'pip'),
        'prokka_gene_fusion':('PE21', 'PE21'),
    }
    ratt = features[source_genome[pair]][pairs[pair][0]]['ratt']
    abinit = features[source_genome[pair]][pairs[pair][1]]['abinit']

    record_sequence = list(SeqIO.parse(f'data/{source_genome[pair]}.fasta', 'fasta'))[0]
    for f in [ratt, abinit]:
        f.references = {record_sequence.id: record_sequence.seq}
        for part in f.location.parts:
            part.ref = record_sequence.id
    annomerge.ref_annotation = ref_features[ref_genome[pair]]
    annomerge.genetic_code = 11

    config.hybran_tmp_dir = tmp_path
    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[pair]}.fasta', 'fasta'))[0].seq
    ref_sequence = SeqIO.read('data/H37Rv.fasta', 'fasta').seq

    expected = {
        'ratt_better': (
            False, True,
            'worse_ref_correspondence',
            "RATT annotation more accurately named and delineated compared to the ab initio annotation.",
        ),
        'ratt_better_coverage': (
            False, True,
            'worse_ref_correspondence',
            "RATT annotation more accurately named and delineated compared to the ab initio annotation.",
        ),
        'pseudo_vs_nonpseudo': (
            False, True,
            'pseudo',
            "Non-pseudo RATT annotation takes precedence over the pseudo ab initio annotation.",
        ),
        'overlapping_unnamed': (
            False, True,
            'unnamed',
            "Unnamed gene and conflicts (overlapping in-frame) with named rival annotation.",
        ),
        'abinit_better': (
            True, False,
            'worse_ref_correspondence',
            "Ab initio annotation more accurately named and delineated compared to the RATT annotation.",
        ),
        'overlapping_different_names_ratt_better': (
            False, True,
            'putative_misannotation',
            "Has no reference-corresponding coordinates, while rival feature has a reference-corresponding start, and both share the same stop position.",
        ),
        'overlapping_different_names_abinit_better': (
            True, False,
            'putative_misannotation',
            "Has no reference-corresponding stop, while rival feature does, and both share the same stop position."
        ),
        'ratt_join_vs_prokka_bad_start': (
            True, False,
            'internal_stop',
            "The ab initio annotation is favored over the RATT annotation because it doesn't contain any internal stops.",
        ),
        'prokka_gene_fusion': (
            False, True,
            'worse_ref_correspondence',
            "RATT annotation more accurately named and delineated compared to the ab initio annotation.",
        ),
    }


    assert annomerge.check_inclusion_criteria(
        ratt_annotation=ratt,
        abinit_annotation=abinit,
    ) == expected[pair]
