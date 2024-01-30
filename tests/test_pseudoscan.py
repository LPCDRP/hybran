import pytest

from hybran import (
    annomerge,
    config,
    demarcate,
    pseudoscan,
)
from hybran.bio import SeqIO
from hybran.config import cnf

from .data_features import *


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
    ['compound_location', 95, 95, True],
])
def test_pseudoscan(feature_type, seq_ident, seq_covg, attempt_rescue, tmp_path):
    ref_genome = defaultdict(lambda :'H37Rv')
    ref_genome['compound_location'] = 'PAO1_107'
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
        'compound_location': 'PAK',
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
        'compound_location': features[source_genome['compound_location']]['PA3701']['ratt'],
    }

    feature = test_features[feature_type]
    config.hybran_tmp_dir = tmp_path
    ref_feature = ref_features[ref_genome[feature_type]][
        annomerge.key_ref_gene(test_features[feature_type].source, test_features[feature_type].qualifiers['gene'][0])
    ]
    record_sequence = list(SeqIO.parse(f'data/{source_genome[feature_type]}.fasta', 'fasta'))[0]
    for part in test_features[feature_type].location.parts:
        part.ref = record_sequence.id
    test_features[feature_type].references = {record_sequence.id: record_sequence.seq}
    cnf.genetic_code = 11
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
        'compound_location': [
            False,
            FeatureLocation(1355921, 1355992, strand=1, ref='refseq|NZ_LR657304.1|chromosome1') \
            + FeatureLocation(1355993, 1357017, strand=1, ref='refseq|NZ_LR657304.1|chromosome1'),
        ],
    }
    results = pseudoscan.pseudoscan(feature, ref_feature, seq_ident, seq_covg, attempt_rescue)
    assert [results, feature.location] == expected[feature_type]

