from collections import defaultdict
from copy import deepcopy

import pytest

from hybran import (
    annomerge,
    config,
    demarcate,
    fusionfisher,
)
from hybran.designator import key_ref_gene
from hybran.util import keydefaultdict

from .data_features import *


@pytest.mark.parametrize('gene_list', [
    'misannotation_false_delayed_stop',
    'partial_fusion_diff_start',
    'redundant_double_partial_fusion',
    'misannotation_both_nonpseudo',
    'misannotation_one_pseudo',
    'idempotence',
])
def test_fusionfisher(gene_list):
    source_genomes = {
        'misannotation_false_delayed_stop':'SEA17020030',
        'partial_fusion_diff_start':'1-0006',
        'redundant_double_partial_fusion':'AZ20',
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
        'partial_fusion_diff_start': [
            source_features['Rv0074']['ratt'],
            source_features['Rv0071']['ratt'],
        ],
        'redundant_double_partial_fusion': [
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
        'redundant_double_partial_fusion': 'nissle-hybrid',
        'misannotation_both_nonpseudo': 'nissle-hybrid',
        'misannotation_one_pseudo': 'nissle-hybrid',
        'idempotence': 'nissle-hybrid',
    })
    source_genome = source_genomes[gene_list]

    record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
    config.cnf.genetic_code = 11
    annomerge.ref_annotation = keydefaultdict(annomerge.ref_fuse)
    annomerge.ref_annotation.update(ref_features[ref_genome[gene_list]])
    for f in inputs[gene_list]:
        f.references = {record_sequence.id: record_sequence.seq}
        for part in f.location.parts:
            part.ref = record_sequence.id
        (f.rcs, f.rce) = demarcate.coord_check(
            f, annomerge.ref_annotation[key_ref_gene(f.source, f.qualifiers['gene'][0])]
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
        'partial_fusion_diff_start': (
            [ ],
            [ ],
            [{
                'feature':source_features['Rv0074']['ratt'],
                'evid':'combined_annotation',
                'remark':"Apparent partial gene fusion.",
            }],
        ),
        'redundant_double_partial_fusion': (
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
    elif gene_list == 'partial_fusion_diff_start':
        fusion_result = fusionfisher.fusion_upgrade(
            base=deepcopy(source_features['Rv0071']['ratt']),
            upstream=source_features['Rv0071']['ratt'],
            downstream=source_features['Rv0074']['ratt'],
            update_location=True
        )
        for outlist in expected['partial_fusion_diff_start'][0:2]:
            outlist.append(fusion_result)

    # set the rival feature as the passing one. for these test cases, we're always looking at pairs
    if expected[gene_list][2]:
        expected[gene_list][2][0]['superior'] = expected[gene_list][0][0]

    assert fusionfisher.fusionfisher(
        deepcopy(inputs[gene_list]),
        annomerge.ref_annotation,
    ) == expected[gene_list]
