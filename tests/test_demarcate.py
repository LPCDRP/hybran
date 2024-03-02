import os
from collections import defaultdict

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pytest

from hybran import (
    annomerge,
    config,
    demarcate,
    extractor,
    pseudoscan,
    designator,
)
from hybran.bio import SeqIO
from hybran.util import keydefaultdict

from .data_features import *

@pytest.mark.parametrize('case', [
    'selenocysteine',
    'no_stop_in_ref',
    'no_stop_in_feature',
    'multi_transl_except',
])
def test_update_transl_except(case, tmp_path):
    config.cnf.genetic_code = 11
    config.hybran_tmp_dir = tmp_path

    if case == 'selenocysteine':
        ref_genome = 'PAO1_107'
        source_genome = 'PAK'
        record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
        feature = deepcopy(features[source_genome]['fdnG']['ratt'])
        feature.location.ref = record_sequence.id
        feature.references = {record_sequence.id: record_sequence.seq}
        ref_feature = deepcopy(ref_features[ref_genome][f"{feature.source}@@@{feature.qualifiers['gene'][0]}"])
        expected = deepcopy(feature)

        expected.qualifiers['transl_except'] = ['(pos:5520855..5520853,aa:Sec)']

    elif case == 'no_stop_in_ref':
        ref_genome = 'PAO1_107'
        source_genome = 'PAK'
        record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
        feature = deepcopy(features[source_genome]['fdnG']['ratt'])
        feature.location.ref = record_sequence.id
        feature.references = {record_sequence.id: record_sequence.seq}

        #Modify the selenocysteine codon in the reference
        ref_feature = deepcopy(ref_features[ref_genome][f"{feature.source}@@@{feature.qualifiers['gene'][0]}"])
        mutant_refseq = list(ref_feature.references[feature.source])

        #mutant_refseq[5401117:5401120] corresponds to '(pos:complement(5401118..5401120),aa:Sec)'
        #which becomes 'TCA' reverse complemented to 'TGA'
        mutant_refseq[5401119] = 'T' #becomes 'AGA'
        mutant_refseq = Seq(''.join(mutant_refseq))
        ref_feature.location.ref = feature.source
        ref_feature.references[feature.source] = mutant_refseq
        expected = deepcopy(feature)

        expected.qualifiers.pop('transl_except')

    elif case == 'no_stop_in_feature':
        ref_genome = 'PAO1_107'
        source_genome = 'PAK'
        record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]

        #Modify the selenocysteine codon in the feature
        mutant_seq = list(record_sequence)

        #mutant_seq[5520853:5520856] corresponds to '(pos:5520855..5520853,aa:Sec)'
        #which becomes 'TCA' reverse complemented to 'TGA'
        mutant_seq[5520855] = 'T' #becomes 'AGA'
        mutant_seq = Seq(''.join(mutant_seq))
        feature = deepcopy(features[source_genome]['fdnG']['ratt'])
        feature.location.ref = record_sequence.id
        feature.references = {record_sequence.id: mutant_seq}
        ref_feature = deepcopy(ref_features[ref_genome][f"{feature.source}@@@{feature.qualifiers['gene'][0]}"])
        expected = deepcopy(feature)

        expected.qualifiers.pop('transl_except')

    elif case == 'multi_transl_except':
        ref_genome = 'PAO1_107'
        source_genome = 'PAK'
        record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]

        #Modify the Pyrrolysine codon in the feature
        #For reference mutant_seq[5520853:5520863] and mutant_refseq[5401117:5401127] is
        #['T', 'C', 'A', 'G', 'A', 'C', 'A', 'C', 'G', 'C']

        mutant_seq = list(record_sequence)
        mutant_seq[5520859:5520862] = 'CTA' #reverse complemented to 'TAG' (originally 'ACG')
        mutant_seq = Seq(''.join(mutant_seq))
        feature = deepcopy(features[source_genome]['fdnG']['ratt'])
        feature.ref = record_sequence.id
        feature.references = {record_sequence.id: mutant_seq}

        #Modify the Pyrrolysine codon in the reference
        ref_feature = deepcopy(ref_features[ref_genome][f"{feature.source}@@@{feature.qualifiers['gene'][0]}"])
        mutant_refseq = list(ref_feature.references[feature.source])
        mutant_refseq[5401123:5401126] = 'CTA' #reverse complemented to 'TAG' (originally 'ACG')
        mutant_refseq = Seq(''.join(mutant_refseq))
        designator.append_qualifier(
            ref_feature.qualifiers,
            'transl_except',
            '(pos:5401124..5401126,aa:Pyl)',
        )

        ref_feature.ref = feature.source
        ref_feature.references[feature.source] = mutant_refseq

        expected = deepcopy(feature)
        expected.qualifiers['transl_except'] = [
            '(pos:5520855..5520853,aa:Sec)',
            '(pos:5520861..5520859,aa:Pyl)'
        ]


    for f in [feature, expected]:
        f.qualifiers['translation'] = [ str(f.translate(
            record_sequence.seq,
            table=config.cnf.genetic_code,
            to_stop=True,
            cds=False,
        )) ]

    pseudoscan.pseudoscan(expected, ref_feature, 95, 95)
    pseudoscan.pseudoscan(feature, ref_feature, 95, 95)
    assert feature == expected

@pytest.mark.parametrize('case,circular', [
    ['softball', False],
    ['minus_strand', False],
    ['compound', False],
])
def test_stopseeker(case, circular):
    if case == 'softball':
        source_genome = '1-0006'
        gene = features[source_genome]['dnaA']['ratt']
        expected = deepcopy(gene)
        # artificially shorten the gene to let stopseeker fix it
        gene.location._end -= 33
    elif case == 'minus_strand':
        source_genome = '1-0006'
        gene = features[source_genome]['Rv1075c']['ratt']
        expected = deepcopy(gene)
        gene.location._start += 8
    elif case == 'compound':
        source_genome = 'PAK'
        gene = features[source_genome]['PA3701']['ratt']
        expected = deepcopy(gene)
        gene.location.parts[-1]._end -= 33


    config.cnf.genetic_code = 11
    record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
    for part in gene.location.parts + expected.location.parts:
        part.ref = record_sequence.id
    gene.references = {record_sequence.id: record_sequence.seq}

    assert demarcate.stopseeker(gene, circular) == expected

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
    test_features[feature_type].location.ref = record_sequence.id
    test_features[feature_type].references = {record_sequence.id: record_sequence.seq}
    config.cnf.genetic_code = 11
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
            'og_de':False,
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
    results = demarcate.coord_check(feature, ref_feature, fix_start=fix_start, fix_stop=fix_stop, seek_stop=seek_stop)
    assert [
        [results, feature.location],
        feature.og.de, feature.corr.de
    ] == [
        expected[feature_type]['results'],
        expected[feature_type]['og_de'],
        expected[feature_type]['corr_de'],
    ]
