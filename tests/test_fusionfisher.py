from collections import defaultdict
from copy import deepcopy

from Bio.SeqFeature import FeatureLocation
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
    'whole_gene_fusion',
    'misannotation_false_delayed_stop',
    'partial_fusion_diff_start',
    'redundant_double_partial_fusion',
    'misannotation_both_nonpseudo',
    'misannotation_one_pseudo',
    'idempotence',
    'consecutive_fusions',
])
def test_fusionfisher(gene_list):

    if gene_list == 'consecutive_fusions':
        pytest.skip("incomplete test case")

    source_genomes = {
        'whole_gene_fusion':'1-0006',
        'misannotation_false_delayed_stop':'SEA17020030',
        'partial_fusion_diff_start':'1-0006',
        'redundant_double_partial_fusion':'AZ20',
        'misannotation_both_nonpseudo': 'AZ20',
        'misannotation_one_pseudo': 'AZ20',
        'idempotence': 'AZ20',
        'consecutive_fusions': 'AZ20_reannotatedref',
    }
    # create a dummy feature dictionary for the test case that is not in focus in the current invocation to avoid a KeyError when loading all expected inputs and results
    source_features = defaultdict(lambda :defaultdict(dict))
    source_features.update(features[source_genomes[gene_list]])
    inputs = {
        'whole_gene_fusion': [
            source_features['Rv0325']['ratt'],
            source_features['Rv0326']['ratt'],
        ],
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
        'consecutive_fusions': [
            source_features['ECOLIN_26716']['ratt'],
            source_features['ECOLIN_26715']['ratt'],
            source_features['ECOLIN_26714']['ratt_corrected'],
            source_features['ECOLIN_26713']['ratt'],
        ]
    }
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome.update({
        'redundant_double_partial_fusion': 'nissle-hybrid',
        'misannotation_both_nonpseudo': 'nissle-hybrid',
        'misannotation_one_pseudo': 'nissle-hybrid',
        'idempotence': 'nissle-hybrid',
        'consecutive_fusions': 'ECOLIN',
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

    # outputs are [outlist, fusions, rejects]
    expected = {
        'whole_gene_fusion': (
            deepcopy(inputs['whole_gene_fusion']),
            [ ], # populated below
            [ ],
        ),
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
        'consecutive_fusions': (
            deepcopy(inputs['consecutive_fusions']),
            [],
            [],
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
        component1 = deepcopy(source_features['Rv0071']['ratt'])
        component1.location = FeatureLocation(81235, 81254, strand=1, ref='1')
        component2 = deepcopy(source_features['Rv0074']['ratt'])
        component2.location = FeatureLocation(81327, 82276, strand=1, ref='1')
        fusion_result.fusion_type = 'partial'
        fusion_result.fusion_components = [
            component1,
            component2,
        ]
        for outlist in expected['partial_fusion_diff_start'][0:2]:
            outlist.append(fusion_result)
    elif gene_list == 'whole_gene_fusion':
        expected['whole_gene_fusion'][0][1].qualifiers['note'] = [
            "Upstream gene Rv0325|Rv0325 conjoins with this one."
        ]
        component1 = deepcopy(expected['whole_gene_fusion'][0][0])
        expected['whole_gene_fusion'][0][0].qualifiers['gene'] = ['Rv0325::Rv0326']
        expected['whole_gene_fusion'][0][0].fusion_type = "whole"
        expected['whole_gene_fusion'][0][0].fusion_components = [
            component1,
            expected['whole_gene_fusion'][0][1],
        ]
        component1.location = FeatureLocation(392626, 392851, strand=1, ref='1')
        expected['whole_gene_fusion'][1].append(expected['whole_gene_fusion'][0][0])
    elif gene_list == 'consecutive_fusions':
        # TODO: check all these fusion component coordinates. see if they make sense
        #
        # NOTE**  - NISSLEORF0154 fusions occur twice. There seem to be two consecutive copies that got fused
        #           but our practice is to not name them twice (...::NISSLEORF0154::NISSLEORF0154)
#         [
# 1[2651439:2658051](-)
# Fusion type: whole
# Components:
# NISSLEORF0154|1[2656502:2658051](-)
# ECOLIN_26716|NISSLEORF0143|1[2651439:2656341](-)
# Qualifiers:
# {'locus_tag': ['ECOLIN_26714'], 'gene': ['NISSLEORF0154::NISSLEORF0143'], 'note': ['Upstream gene ECOLIN_26714|NISSLEORF0154 conjoins with this one.']}
# ,
# 1[2651439:2659227](-)
# Fusion type: whole
# Components:
# NISSLEORF0154|1[2656502:2659227](-)
# ECOLIN_26714|NISSLEORF0154::NISSLEORF0143|1[2651439:2658051](-)
# Qualifiers:
# {'locus_tag': ['ECOLIN_26714'], 'gene': ['NISSLEORF0154::NISSLEORF0143'], 'pseudo': [''], 'note': ['Hybran/Pseudoscan:description:Locus does not have reference-corresponding start | No internal stop codons and ends with a valid stop codon | Locus has valid reading frame | Locus has a valid alternative start site', 'Hybran/Pseudoscan:evidence:no_rcc', 'Hybran/Pseudoscan:barcode:D31;VS1;VE1;RCS0;RCE1;BOK.', 'Upstream gene ECOLIN_26713|NISSLEORF0188 conjoins with this one.']}
# ,
# 1[2651439:2659977](-)
# Fusion type: whole
# Components:
# NISSLEORF0188|1[2659130:2659977](-)
# ECOLIN_26714|NISSLEORF0154::NISSLEORF0143|1[2651439:2659227](-)
# Qualifiers:
# {'locus_tag': ['ECOLIN_26713'], 'gene': ['NISSLEORF0188::NISSLEORF0154::NISSLEORF0143']}
# ]
        expected['consecutive_fusions'][1] = expected['consecutive_fusions'][0][1:]
        # first fusion
        expected['consecutive_fusions'][1][0].qualifiers['note'] = [
            'Upstream gene ECOLIN_26714|NISSLEORF0154 conjoins with this one.',
        ]
        expected['consecutive_fusions'][1][0].fusion_type = "whole"
        componentX = deepcopy(expected['consecutive_fusions'][1][0])
        (componentX.location._start, componentX.location._end)  = (2656502, 2658051)
        expected['consecutive_fusions'][1][0].components = [
            componentX,
            expected['consecutive_fusions'][0][0],
        ]
        # second fusion
        expected['consecutive_fusions'][1][1].qualifiers['note'] = [
            'Upstream gene ECOLIN_26713|NISSLEORF0188 conjoins with this one.',
        ]
        expected['consecutive_fusions'][1][1].fusion_type = "whole"
        componentY = deepcopy(expected['consecutive_fusions'][1][1])
        (componentY.location._start, componentY.location._end)  = (2656502, 2658051)
        expected['consecutive_fusions'][1][1].components = [
            expected['consecutive_fusions'][0][0],
            #....
        ]
        # third fusion

    # set the rival feature as the passing one. for these test cases, we're always looking at pairs
    if expected[gene_list][2]:
        expected[gene_list][2][0]['superior'] = expected[gene_list][0][0]

    assert fusionfisher.fusionfisher(
        deepcopy(inputs[gene_list]),
        annomerge.ref_annotation,
    ) == expected[gene_list]
