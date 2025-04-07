from collections import defaultdict
import os

from Bio.SeqFeature import (
    SeqFeature,
    SimpleLocation,
    FeatureLocation,
    ExactPosition,
    CompoundLocation,
)
import pytest

from hybran import (
    annomerge,
    config,
)
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
    config.cnf.genetic_code = 11

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
