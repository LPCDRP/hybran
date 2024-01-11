from collections import defaultdict, OrderedDict
from copy import deepcopy
import os

from Bio.SeqFeature import SeqFeature, SimpleLocation, FeatureLocation, ExactPosition, CompoundLocation

import pytest

from hybran import annomerge
from hybran import config
from hybran import demarcate
from hybran import pseudoscan
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
    annomerge.genetic_code = demarcate.genetic_code = pseudoscan.genetic_code = 11
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
    'hybrid_fusion_diff_start',
    'redundant_double_hybrid_fusion',
    'misannotation_both_nonpseudo',
    'misannotation_one_pseudo',
    'idempotence',
])
def test_fusionfisher(gene_list):
    source_genomes = {
        'misannotation_false_delayed_stop':'SEA17020030',
        'hybrid_fusion_diff_start':'1-0006',
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
            source_features['Rv3382c']['ratt'],
            source_features['Rv3383c']['ratt'],
        ],
        'hybrid_fusion_diff_start': [
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
    annomerge.genetic_code = demarcate.genetic_code = 11
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[gene_list]])
    for f in inputs[gene_list]:
        f.references = {record_sequence.id: record_sequence.seq}
        for part in f.location.parts:
            part.ref = record_sequence.id
        (f.rcs, f.rce) = annomerge.coord_check(
            f, annomerge.ref_annotation[annomerge.key_ref_gene(f.source, f.qualifiers['gene'][0])]
        )

    expected = {
        'misannotation_false_delayed_stop': (
            [ source_features['Rv3383c']['ratt'] ],
            [ ],
            [{
                'feature':source_features['Rv3382c']['ratt'],
                'evid':'putative_misannotation',
                'remark':'Has no reference-corresponding stop, while rival feature does, and both share the same stop position.',
            }],
        ),
        'hybrid_fusion_diff_start': (
            [ ],
            [ ],
            [{
                'feature':source_features['Rv0074']['ratt'],
                'evid':'combined_annotation',
                'remark':"Apparent hybrid fusion gene.",
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
    elif gene_list == 'hybrid_fusion_diff_start':
        fusion_result = annomerge.fusion_upgrade(
            base=deepcopy(source_features['Rv0071']['ratt']),
            upstream=source_features['Rv0071']['ratt'],
            downstream=source_features['Rv0074']['ratt'],
            update_location=True
        )
        for outlist in expected['hybrid_fusion_diff_start'][0:2]:
            outlist.append(fusion_result)

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
    annomerge.genetic_code = demarcate.genetic_code = 11

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
            'forfeit',
            "Equally valid call. RATT annotation is favored due to synteny.",
        ),
        'overlapping_different_names_abinit_better': (
            True, False,
            'worse_ref_correspondence',
            "Ab initio annotation more accurately named and delineated compared to the RATT annotation.",
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
