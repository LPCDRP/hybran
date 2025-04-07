from collections import defaultdict
from copy import deepcopy
import os

from Bio.SeqFeature import (
    SeqFeature,
    SimpleLocation,
    FeatureLocation,
    ExactPosition,
    CompoundLocation,
)
import pytest

from hybran import ratt

from .data_features import *


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
    Rv2434c =  SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(2728810), ExactPosition(2729608), strand=-1), FeatureLocation(ExactPosition(2728376), ExactPosition(2728805), strand=-1)], 'join'), type='CDS')
    Rv2434c.qualifiers = dict(locus_tag=["Rv2434c"], codon_start=['1'], transl_table='11')
    # Rv0739 - compound CDS with internal stop codons in this isolate
    Rv0739 = SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(829165), ExactPosition(829201), strand=1), FeatureLocation(ExactPosition(829204), ExactPosition(829993), strand=1)], 'join'), type='CDS')
    Rv0739.qualifiers = dict(locus_tag=["Rv0739"], codon_start=['1'], transl_table='11')
    Rv3020c = SeqFeature(FeatureLocation(ExactPosition(3371022), ExactPosition(3371217), strand=-1), type='CDS')
    Rv3020c.qualifiers = dict(locus_tag=["Rv3020c"], codon_start=['1'], transl_table='11')

    ref_annotation = ref_features[ref_genome[case]]
    ratt.genetic_code = 11

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

    record_sequence = SeqIO.read("data/1-0009.fasta",format="fasta").seq
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

    assert ratt.isolate_valid_ratt_annotations(
        feature_list,
        ref_annotation,
        reference_locus_list,
        seq_ident,
        seq_covg,
        ratt_enforce_thresholds=filter,
    ) == expected

#                                             'pseudo' in feature_list[3].qualifiers.keys() and feature_list[3].qualifiers['pseudo'] == [''] # Rv2434c
