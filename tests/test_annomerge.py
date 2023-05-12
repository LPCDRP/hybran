from collections import defaultdict, OrderedDict
from copy import deepcopy
import os

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, FeatureLocation, ExactPosition, CompoundLocation

import pytest

from hybran import annomerge
from hybran import config

from .data_features import *



def test_ref_fuse():
    annomerge.ref_annotation = annomerge.keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update({'@@@'.join(['H37Rv',k]):v for k,v in ref_features['H37Rv'].items()})

    assert annomerge.ref_annotation['H37Rv::H37Rv@@@PE_PGRS50::PE_PGRS49'].location == CompoundLocation([
        SimpleLocation(ExactPosition(3738157), ExactPosition(3742774), strand=-1),
        SimpleLocation(ExactPosition(3736983), ExactPosition(3738000), strand=-1)
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
    }

    assert annomerge.overlap_inframe(pairs[pair][0], pairs[pair][1]) == expected[pair]

@pytest.mark.parametrize('gene_list', [
    'complementary_fragments',
    'complementary_fragments_one_unnamed',
    'fails_final_coord_check_inframe_overlap',
#    'seemingly_complete_fragment',
#    'independent_fragments',
])
@pytest.mark.skipif(not os.path.isfile("data/1-0006.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/2-0031.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.fasta"), reason="test reference genome sequence not available")
def test_fissionfuser(gene_list):
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
        'independent_fragments': [
        ],
    }
    source_genomes = {
        'complementary_fragments':'1-0006',
        'complementary_fragments_one_unnamed':'1-0006',
        'fails_final_coord_check_inframe_overlap':'2-0031',
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    source_genome = source_genomes[gene_list]
    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0].seq
    annomerge.genetic_code = 11
    annomerge.ref_annotation = ref_features[ref_genome[gene_list]]
    annomerge.corrected_orf_report = []

    expected = {
        'complementary_fragments': ([
            SeqFeature(
                FeatureLocation(ExactPosition(2275540), ExactPosition(2277261), strand=-1), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02174'],
                    'pseudo': [''],
                }
            )],
            [(inputs['complementary_fragments'][0], 'L_02173:dosT combined with L_02174:dosT: overlapping_inframe')]
        ),
        'complementary_fragments_one_unnamed': ([
            SeqFeature(
                FeatureLocation(ExactPosition(1104310), ExactPosition(1105056), strand=1), type='CDS', qualifiers={
                    'gene': ['Rv0986'],
                    'locus_tag': ['L_01054'],
                    'pseudo': [''],
                }
            )],
            [(inputs['complementary_fragments_one_unnamed'][0], 'L_01053:L_01053 combined with L_01054:Rv0986: complementary_fragments')]
        ),
        'fails_final_coord_check_inframe_overlap': ([
            SeqFeature(
                FeatureLocation(ExactPosition(3522323), ExactPosition(3523418), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03352'],
                    'pseudo': [''],
                }
            )],
            [(inputs['fails_final_coord_check_inframe_overlap'][0], 'L_03351:Rv3327 combined with L_03352:Rv3327: overlapping_inframe')]
        ),
    }

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
])
def test_fusionfisher(gene_list):
    source_genomes = {
        'misannotation_false_delayed_stop':'1-0006',
        'redundant_double_hybrid_fusion':'AZ20',
        'misannotation_both_nonpseudo': 'AZ20',
        'misannotation_one_pseudo': 'AZ20',
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
    }
    expected = {
        'misannotation_false_delayed_stop': (
            [ source_features['Rv0074']['ratt'] ],
            [ ],
            [ (source_features['Rv0071']['ratt'], "putative misannotation: has no reference-corresponding stop, while Rv0074:Rv0074 does, and both share the same stop position.") ]
        ),
        'redundant_double_hybrid_fusion': (
            [ source_features['AZ20_03933']['ratt'] ],
            [ ],
            [ (source_features['AZ20_03933']['prokka'], "Redundant annotation with ECOLIN_01320:ORF0033::ECOLIN_01320") ],
        ),
        'misannotation_both_nonpseudo': (
            [ source_features['ECOLIN_18975']['ratt'] ],
            [ ],
            [ (source_features['ECOLIN_18965']['ratt'] , "Both ECOLIN_18975:ECOLIN_18975 and ECOLIN_18965:ECOLIN_18965 have reference-corresponding start codons and have reference-corresponding stop codons. Longer feature ECOLIN_18975:ECOLIN_18975 favored.") ],
        ),
        'misannotation_one_pseudo': (
            [source_features['ECOLIN_25305']['ratt']],
            [ ],
            [(source_features['garD']['ratt'], "putative misannotation: has no reference-corresponding coordinates, while ECOLIN_25305:ECOLIN_25305 has a reference-corresponding start, and  both share the same stop position.")],
        ),
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome.update({
        'redundant_double_hybrid_fusion': 'nissle-hybrid',
        'misannotation_both_nonpseudo': 'nissle-hybrid',
        'misannotation_one_pseudo': 'nissle-hybrid',
    })
    source_genome = source_genomes[gene_list]

    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0].seq
    annomerge.genetic_code = 11
    annomerge.ref_annotation = annomerge.keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[gene_list]])

    assert annomerge.fusionfisher(
        inputs[gene_list],
    ) == expected[gene_list]

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
    annomerge.liftover_annotation(abinit, ref, False, inference='alignment:blastp')
    assert abinit.qualifiers == expected

@pytest.mark.parametrize('feature_type,fix_start,fix_stop', [
    ['abinit_start_bad_minus', True, False],
    ['bad_translation', True, False],
    ['tricky_found_low', True, False],
    ['good_start_stop_deletion', True, False],
    ['bad_mismatch_check', True, True],
    ['bad_mismatch_check2', True, True],
    ['ref_start_frameshift', True, False],
    ['bad_start_good_padding', True, False],
    ['ratt_pseudo_pgrs', True, True],
    ['same_start_alt_stop_1', False, False],
    ['same_start_alt_stop_2', False, False],
    ['same_start_alt_stop_2_fix', True, True],
    ['bad_start_stop_nofix_pseudo', False, False],
    ['bad_start_stop_fix_pseudo', True, False],
    ['inverted_join_ecoli', True, True],
    ['gene_fusion', True, True],
    ['end_greater_than_start', True, True],
])
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.gbk"), reason="test reference annotation not available")
@pytest.mark.skipif(not os.path.isfile("data/nissle-hybrid.gbk"), reason="test reference annotation not available")
def test_coord_check(feature_type, fix_start, fix_stop):
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
        'inverted_join_ecoli': 'AZ20',
        'gene_fusion': '1-0006',
        'end_greater_than_start': 'SEA08151',

    }
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome['inverted_join_ecoli'] = 'nissle-hybrid'

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
    }

    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[feature_type]}.fasta', 'fasta'))[0].seq
    annomerge.genetic_code = 11
    annomerge.corrected_orf_report = []
    annomerge.ref_annotation = annomerge.keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[feature_type]])

    feature = test_features[feature_type]
    ref_feature = annomerge.ref_annotation[
        annomerge.key_ref_gene(test_features[feature_type].source, test_features[feature_type].qualifiers['gene'][0])
    ]

    expected = {
        'abinit_start_bad_minus': [(True, True), FeatureLocation(3548089, 3548542, strand=-1)],
        'bad_translation':[(True, True), FeatureLocation(4230767, 4231754, strand=1)],
        'tricky_found_low':[(False, True), FeatureLocation(2636045, 2637140, strand=-1)],
        'good_start_stop_deletion':[(True, True), FeatureLocation(4239393, 4240378, strand=1)],
        'bad_mismatch_check':[(True, False), FeatureLocation(2112334, 2114834, strand=1)],
        'bad_mismatch_check2':[(False, False), FeatureLocation(3707086, 3709176, strand =1)],
        'ref_start_frameshift':[(True, True), FeatureLocation(3374234, 3375312, strand=-1)],
        'bad_start_good_padding':[(False, True), FeatureLocation(1927378, 1927894, strand=-1)],
        'ratt_pseudo_pgrs':[(True,False), FeatureLocation(3741108, 3746955, strand=-1)],
        'same_start_alt_stop_1':[(False, True), FeatureLocation(3182302, 3183397, strand=-1)],
        'same_start_alt_stop_2':[(True, False), FeatureLocation(3182302, 3183397, strand=-1)],
        'same_start_alt_stop_2_fix':[(True, True), FeatureLocation(3182570, 3183397, strand=-1)],
        'bad_start_stop_nofix_pseudo': [(False, False), FeatureLocation(1217413, 1217872, strand=1)],
        'bad_start_stop_fix_pseudo': [(True, False), FeatureLocation(1217428, 1217872, strand=1)],
        'inverted_join_ecoli': [(False, False), FeatureLocation(3714209, 3716770, strand=-1)],
        'gene_fusion': [(True, True), FeatureLocation(3741108, 3746955, strand=-1)],
        'end_greater_than_start': [(False, False), FeatureLocation(4064133, 4064322, strand=1)],
    }
    results = annomerge.coord_check(feature, ref_feature, fix_start=fix_start, fix_stop=fix_stop)
    assert [results, feature.location] == expected[feature_type]

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
    }

    feature = test_features[feature_type]
    config.hybran_tmp_dir = tmp_path
    ref_feature = ref_features[ref_genome[feature_type]][
        annomerge.key_ref_gene(test_features[feature_type].source, test_features[feature_type].qualifiers['gene'][0])
    ]
    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[feature_type]}.fasta', 'fasta'))[0].seq
    annomerge.genetic_code = 11
    annomerge.corrected_orf_report = []
    expected = {
        'small_badstop_fix_pseudo': [True, FeatureLocation(66693, 67032, strand=-1)],
        'small_badstart_fix_nopseudo': [False, FeatureLocation(712762, 713308, strand=1)],
        'small_badstop_fix_pseudo_frameshift': [True, FeatureLocation(1202098, 1203042, strand=-1)],
        'frameshift_alt_stop_pseudo' : [True, FeatureLocation(192566, 193259, strand=-1)],
        'good_start_stop_frameshift_pseudo' : [True, FeatureLocation(366753, 376299, strand=-1)],
        'bad_start_stop_nofix_pseudo': [True, FeatureLocation(1217428, 1217872, strand=1)],
        'fix_stop_valid_broken_stop': [True, FeatureLocation(2339465, 2341436, strand=1)],
        'sensitive_padding_fix_pseudo': [True, FeatureLocation(2342175, 2342617, strand=-1)],
        'deletion_in_middle_fix_pseudo': [True, FeatureLocation(2397532, 2398106, strand=-1)],
        'fix_start_stop_nonpseudo': [False, FeatureLocation(2552370, 2553750, strand=1)],
        'good_start_stop_fix_pseudo': [True, FeatureLocation(2735891, 2736310, strand=1)],
        'inframe_deletion_in_middle': [False, FeatureLocation(3730264, 3741334, strand=-1)],
        'good_blast_still_repairable': [False, FeatureLocation(0, 1524, strand=1)],
        'start_correction_induces_delayed_stop': [True, FeatureLocation(1370377, 1371385, strand=-1)],
        'start_correction_induces_delayed_stop2': [True, FeatureLocation(2881559, 2882297, strand=1)],
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

@pytest.mark.parametrize("case",[
    'long_unbroken_pseudo',
])
@pytest.mark.skipif(not os.path.isfile("data/1-0009.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.gbk"), reason="test reference annotation not available")
def test_isolate_valid_ratt_annotations(case):
    source_genome = {
        'long_unbroken_pseudo':'1-0006',
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    cases = {
        'long_unbroken_pseudo': 'PE_PGRS50',
    }
    feature_list = [ features[source_genome[case]][cases[case]]['ratt_raw'] ]

    Rv0001 = SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1524), strand=1), type='CDS')
    Rv0001.qualifiers = dict(locus_tag=["Rv0001"], codon_start=['1'], transl_table='11')
    Rv0071 = SeqFeature(FeatureLocation(ExactPosition(81037), ExactPosition(82096), strand=1), type='CDS')
    Rv0071.qualifiers = dict(locus_tag=["Rv0071"], codon_start=['1'], transl_table='11')
    Rv2434c =  SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(2728810), ExactPosition(2729608), strand=-1), FeatureLocation(ExactPosition(2728376), ExactPosition(2728805), strand=-1)], 'join'), type='CDS', location_operator='join')
    Rv2434c.qualifiers = dict(locus_tag=["Rv2434c"], codon_start=['1'], transl_table='11')
    # Rv0739 - compound CDS with internal stop codons in this isolate
    Rv0739 = SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(829165), ExactPosition(829201), strand=1), FeatureLocation(ExactPosition(829204), ExactPosition(829993), strand=1)], 'join'), type='CDS', location_operator='join')
    Rv0739.qualifiers = dict(locus_tag=["Rv0739"], codon_start=['1'], transl_table='11')
    Rv3020c = SeqFeature(FeatureLocation(ExactPosition(3371022), ExactPosition(3371217), strand=-1), type='CDS')
    Rv3020c.qualifiers = dict(locus_tag=["Rv3020c"], codon_start=['1'], transl_table='11')

    annomerge.ref_annotation = ref_features[ref_genome[case]]
    annomerge.genetic_code = 11

    ref_temp_fasta_dict = dict(
        Rv0001  = 'gene-seqs/Rv0001.fasta',
        Rv0071  = 'gene-seqs/Rv0071.fasta',
        Rv2434c = 'gene-seqs/Rv2434c.fasta',
        Rv3020c = 'gene-seqs/Rv3020c.fasta'
        )

    reference_locus_list = ['Rv0001','Rv0071','Rv2434c','Rv3020c']
    filter = False
    if filter:
        seq_ident = seq_covg = 95
    else:
        seq_ident = seq_covg = 0

    annomerge.record_sequence = SeqIO.read("/grp/valafar/data/genomes/1-0009.fasta",format="fasta").seq
    if case == 'long_unbroken_pseudo':
        out_list = deepcopy(feature_list)
        out_list[0].qualifiers['pseudo'] = ['']
        expected = (out_list, {}, [])

    last_expected = {
        True : (
            [Rv2434c, Rv0001],
            {
                "Rv0001": {
                    'iden': 99.803,
                    'scov': 100.0,
                    'qcov': 100.0
                },
            },
            [
                (Rv0739, "Multiple internal stop codons in compound CDS feature."),
                (Rv0071, "No blastp hit to corresponding reference CDS at specified thresholds."),
                (Rv3020c, "No blastp hit to corresponding reference CDS at specified thresholds."),
            ]
        ),
        False : (
            [Rv2434c, Rv0001, Rv0071, Rv3020c],
            {
                "Rv0001": {
                    'iden': 99.803,
                    'scov': 100.0,
                    'qcov': 100.0,
                },
                "Rv0071": {
                    'iden': 45.0,
                    'scov': 7.659574468085106,
                    'qcov': 5.681818181818182,
                },
                "Rv3020c": {
                    'iden': 0.0,
                    'scov': 0.0,
                    'qcov': 0.0,
                },
            },
            [
                (Rv0739, "Multiple internal stop codons in compound CDS feature."),
            ]
        )
    }

    assert annomerge.isolate_valid_ratt_annotations(feature_list,
                                             ref_temp_fasta_dict,
                                             reference_locus_list,
                                             seq_ident, seq_covg) == expected
#                                             'pseudo' in feature_list[3].qualifiers.keys() and feature_list[3].qualifiers['pseudo'] == [''] # Rv2434c


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

    # 1-0006. abinit overlapping fused reference
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
        'overlapping_fused_ref': []
    }
    expected_rejects = {
        'double_overlap': [],
        'overlapping_fused_ref': [(abinit_ppe6, 'duplicate of Rv0305c')]
    }
    expected_overlaps = {
        'double_overlap': defaultdict(list, {loc_triplet(abinit_nissleorf0212): [loc_triplet(nissleorf0212)]}),
        'overlapping_fused_ref': defaultdict(list, {})
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

    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[pair]}.fasta', 'fasta'))[0].seq
    annomerge.ref_annotation = ref_features[ref_genome[pair]]
    annomerge.genetic_code = 11

    config.hybran_tmp_dir = tmp_path
    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[pair]}.fasta', 'fasta'))[0].seq
    ref_sequence = SeqIO.read('data/H37Rv.fasta', 'fasta').seq
    annomerge.ref_prom_fp_dict = annomerge.get_nuc_seq_for_gene(
        [ref_features['H37Rv']['@@@'.join(['H37Rv.NC_000962.3', pairs[pair][0]])]],
        ref_sequence,
    )[0]

    expected = {
        'ratt_better': (
            False, True,
            "RATT annotation Rv0001:dnaA more accurately named and delineated.",
        ),
        'ratt_better_coverage': (
            False, True,
            "RATT annotation Rv1453:Rv1453 more accurately named and delineated.",
        ),
        'pseudo_vs_nonpseudo': (
            False, True,
            "Non-pseudo RATT annotation takes precedence.",
        ),
        'overlapping_unnamed': (False, True, "Hypothetical gene and conflicts (overlapping in-frame) with RATT's Rv0001:dnaA."),
        'abinit_better': (True, False, 'Ab initio annotation L_02383:Rv1718 more accurately named and delineated.'),
        'overlapping_different_names_ratt_better': (
            False, True,
            "putative misannotation: has no reference-corresponding coordinates, while Rv1945:Rv1945 has a reference-corresponding start, and  both share the same stop position.",
        ),
        'overlapping_different_names_abinit_better': (
            True, False,
            "putative misannotation: has no reference-corresponding stop, while L_02335:ORF0004 does, and both share the same stop position."
        ),
        'ratt_join_vs_prokka_bad_start': (
            False, True,
            "Equally valid call, but the more complete RATT annotation is favored."
        ),
        'prokka_gene_fusion': (
            True, False,
            "The ab initio annotation is favored due to having a valid delayed stop."
        ),
    }


    assert annomerge.check_inclusion_criteria(
        ratt_annotation=ratt,
        abinit_annotation=abinit,
    ) == expected[pair]
