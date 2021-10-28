from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, CompoundLocation

import pytest

from hybran import annomerge

@pytest.mark.parametrize("filter",[True, False])
def test_isolate_valid_ratt_annotations(filter):
    Rv0001 = SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1524), strand=1), type='CDS')
    Rv0001.qualifiers = dict(locus_tag=["Rv0001"], codon_start=['1'], transl_table='11')
    Rv0071 = SeqFeature(FeatureLocation(ExactPosition(81037), ExactPosition(82096), strand=1), type='CDS')
    Rv0071.qualifiers = dict(locus_tag=["Rv0071"], codon_start=['1'], transl_table='11')
    Rv2434c =  SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(2728810), ExactPosition(2729608), strand=-1), FeatureLocation(ExactPosition(2728376), ExactPosition(2728805), strand=-1)], 'join'), type='CDS', location_operator='join')
    Rv2434c.qualifiers = dict(locus_tag=["Rv2434c"], codon_start=['1'], transl_table='11')
    Rv3020c = SeqFeature(FeatureLocation(ExactPosition(3371022), ExactPosition(3371217), strand=-1), type='CDS')
    Rv3020c.qualifiers = dict(locus_tag=["Rv3020c"], codon_start=['1'], transl_table='11')
    feature_list = [Rv0001, Rv0071, Rv2434c, Rv3020c]
                   
    
    ref_temp_fasta_dict = dict(
        Rv0001  = 'gene-seqs/Rv0001.fasta',
        Rv0071  = 'gene-seqs/Rv0071.fasta',
        Rv2434c = 'gene-seqs/Rv2434c.fasta',
        Rv3020c = 'gene-seqs/Rv3020c.fasta'
        )

    reference_locus_list = ['Rv0001','Rv0071','Rv2434c','Rv3020c']
    if filter:
        seq_ident = seq_covg = 95
    else:
        seq_ident = seq_covg = 0

    annomerge.record_sequence = SeqIO.read("/grp/valafar/data/genomes/1-0009.fasta",format="fasta").seq
    expected = {True : (
        [Rv0001],
        {"Rv0001": [99.803, 100.0, 100.0]},
        [(Rv2434c, "compound location"),
         (Rv0071, "No blastp hit to corresponding reference CDS at specified thresholds."),
         (Rv3020c, "No blastp hit to corresponding reference CDS at specified thresholds.")]
    ),
                False : ([Rv0001, Rv0071, Rv3020c],
                         {"Rv0001": [99.803, 100.0, 100.0],
                          "Rv0071": [35.714, 5.957446808510639, 3.977272727272727],
                          "Rv3020c": [0.0, 0.0, 0.0]},
                         [(Rv2434c, "compound location")])
                }

    assert annomerge.isolate_valid_ratt_annotations(feature_list,
                                             ref_temp_fasta_dict,
                                             reference_locus_list,
                                             seq_ident, seq_covg) == \
        expected[filter]
