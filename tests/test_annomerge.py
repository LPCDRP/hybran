from collections import defaultdict, OrderedDict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, CompoundLocation

import pytest

from hybran import annomerge
from hybran import config

from .data_features import *


def test_process_split_genes():
    interrupted_by_another_gene = [
        SeqFeature(FeatureLocation(ExactPosition(936376), ExactPosition(937165), strand=-1), type='CDS',
                   qualifiers=OrderedDict([('gene', ['pip']),
                                           ('locus_tag', ['Rv0840c']),
                                           ('note',
                                            ['FunctionalCategory: intermediary metabolism and respiration',
                                             'my_bogus_second_note']),
                                           ('codon_start', ['1']),
                                           ('transl_table', ['11']),
                                           ('product',
                                            ['Probable proline iminopeptidase Pip (prolyl aminopeptidase) '
                                             '(pap)']),
                                           ('inference', ['alignment:RATT']),
                                           ('pseudo', ['']),
                                           ('translation',
                                            ['PRISATIRTARSRSSGGYLLDEPPDMTPSFPRTGVSGHAGAVHVDLWTMNRSVAEMATVAEALALTRFHIFSHSWGGMLAQQYVLDKAPDAVSLTIANSTASIPEFSASLVSLKSCLDVATRSAIDRHEAAGTTHSAEYQAAIRTWNETYLCRTRPWPRELTEAFANMGTEIFETMFGPSDFRIVGNVRDWDVVDRLADIAVPTLLVVGRFDECSPEHMREMQGRIAGSRLEFFESSSHMPFIEEPARFDRVMREFLRLHDI'])])
                   ),
        SeqFeature(FeatureLocation(ExactPosition(937088), ExactPosition(937415), strand=1), type='CDS',
                   qualifiers=OrderedDict([('gene', ['ORF0009']),
                                           ('locus_tag', ['L_00897']),
                                           ('inference',
                                            ['ab initio prediction:Prodigal:002006',
                                             'similar to AA sequence:ref_proteome.fasta:Rv3474:ORF0009']),
                                           ('codon_start', ['1']),
                                           ('transl_table', ['11']),
                                           ('product',
                                            ['putative transposase for insertion element IS6110 (fragment)']),
                                           ('protein_id', ['C:L_00897']),
                                           ('translation',
                                            ['MSGGSSRRYPPELRERAVRMVAEIRGQHDSEWAAISEVARLLGVGCAETVRKWVRQAQVDAGARPGTTTEESAELKRLRRDNAELRRANAILKTASAFFAAELDRPAR']),
                                           ('note', ['This annotation is added from Prokka reference run'])])
                   ),
        SeqFeature(FeatureLocation(ExactPosition(937168), ExactPosition(937360), strand=-1), type='CDS',
                   qualifiers=OrderedDict([('gene', ['pip']),
                                           ('locus_tag', ['Rv0840c']),
                                           ('note',
                                            ['FunctionalCategory: intermediary metabolism and respiration']),
                                           ('codon_start', ['1']),
                                           ('transl_table', ['11']),
                                           ('product',
                                            ['Probable proline iminopeptidase Pip (prolyl aminopeptidase) '
                                             '(pap)']),
                                           ('inference', ['alignment:RATT']),
                                           ('pseudo', ['']),
                                           ('translation',
                                            ['IAFALRNSALSRRKRFSSADSSVVVPGRAPASTCAWRTHLRTVSAQPTPSRRATSLIAAHSESC'])])
                   ),
    ]
    result = annomerge.process_split_genes(interrupted_by_another_gene)
    received = dict(
        locations = [r.location for r in result],
        qualifiers = [r.qualifiers for r in result]
    )
    expected = dict(
        locations = [
            FeatureLocation(ExactPosition(936376), ExactPosition(937360), strand=-1),
            FeatureLocation(ExactPosition(937088), ExactPosition(937415), strand=1)
        ],
        qualifiers = [
            OrderedDict([('gene', ['pip']),
                                           ('locus_tag', ['Rv0840c']),
                                           ('note',
                                            ['FunctionalCategory: intermediary metabolism and respiration',
                                             'my_bogus_second_note']),
                                           ('codon_start', ['1']),
                                           ('transl_table', ['11']),
                                           ('product',
                                            ['Probable proline iminopeptidase Pip (prolyl aminopeptidase) '
                                             '(pap)']),
                                           ('inference', ['alignment:RATT']),
                                           ('pseudo', ['']),
                         ]),
            OrderedDict([('gene', ['ORF0009']),
                                           ('locus_tag', ['L_00897']),
                                           ('inference',
                                            ['ab initio prediction:Prodigal:002006',
                                             'similar to AA sequence:ref_proteome.fasta:Rv3474:ORF0009']),
                                           ('codon_start', ['1']),
                                           ('transl_table', ['11']),
                                           ('product',
                                            ['putative transposase for insertion element IS6110 (fragment)']),
                                           ('protein_id', ['C:L_00897']),
                                           ('translation',
                                            ['MSGGSSRRYPPELRERAVRMVAEIRGQHDSEWAAISEVARLLGVGCAETVRKWVRQAQVDAGARPGTTTEESAELKRLRRDNAELRRANAILKTASAFFAAELDRPAR']),
                                           ('note', ['This annotation is added from Prokka reference run'])
                         ]),
        ]
    )
    assert (
        received['locations'] == expected['locations']
        and received['qualifiers'] == expected['qualifiers']
    )


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

@pytest.mark.parametrize("filter",[True, False])
def test_isolate_valid_ratt_annotations(filter):
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
    feature_list = [Rv0001, Rv0071, Rv0739, Rv2434c, Rv3020c]


    annomerge.genetic_code = 11

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
    expected = {
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
                                             seq_ident, seq_covg) == expected[filter] \
                                             and \
                                             'pseudo' in feature_list[3].qualifiers.keys() and feature_list[3].qualifiers['pseudo'] == [''] # Rv2434c


@pytest.mark.parametrize('pair', [
    'ratt_better',
    'different',
    'abinit_better',
])
def test_check_inclusion_criteria(pair, tmp_path):
    source_genome = {
        'ratt_better':'1-0006',
        'different':'1-0006',
        'abinit_better':'4-0041',
        'corresponding_non_cds':'4-0041',
    }
    pairs = {
        'ratt_better': ('dnaA', 'dnaA'),
        'different': ('dnaA', 'gyrB'),
        'corresponding_non_cds': ('rrf', 'rrf'),
        'similar': (),
        'abinit_better': ('Rv1718', 'Rv1718'),
    }
    ratt = features[source_genome[pair]][pairs[pair][0]]['ratt']
    abinit = features[source_genome[pair]][pairs[pair][1]]['abinit']

    reference_gene_locus_dict = dict(
        dnaA='Rv0001',
        Rv0205='Rv0205',
        rplB='Rv0704',
        Rv1718='Rv1718',
        mamB='Rv2024c',
    )
    reference_locus_gene_dict = dict(
        Rv0001='dnaA',
        Rv0205='Rv0205',
        Rv0704='rplB',
        Rv1718='Rv1718',
        Rv2024c='mamB',
    )

    config.hybran_tmp_dir = tmp_path
    annomerge.record_sequence = list(SeqIO.parse(f'data/{source_genome[pair]}.fasta', 'fasta'))[0].seq
    ref_sequence = SeqIO.read('data/H37Rv.fasta', 'fasta').seq
    annomerge.ref_prom_fp_dict = annomerge.get_prom_for_gene(
        [ref_features['H37Rv'][pairs[pair][0]]],
        ref_sequence,
    )

    expected = {
        'ratt_better': (
            False, True,
            "sequence upstream of RATT's Rv0001 is not mutated.",
        ),
        'different': (True, True, ''),
        'abinit_better': (True, False, 'Ab initio feature L_02383 has better alignment coverage with the reference.'),
    }


    assert annomerge.check_inclusion_criteria(
        ratt_annotation=ratt,
        abinit_annotation=abinit,
        reference_gene_locus_dict=reference_gene_locus_dict,
        reference_locus_gene_dict=reference_locus_gene_dict,
        abinit_blast_results=abinit_blast_results[source_genome[pair]],
        ratt_blast_results=ratt_blast_results[source_genome[pair]],
    ) == expected[pair]
