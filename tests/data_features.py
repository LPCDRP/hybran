from collections import defaultdict
from copy import deepcopy
import os
import pathlib

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, CompoundLocation


ref_features = {}
for ref in [
        'H37Rv',
        'nissle-hybrid',
]:
    ref_file = os.path.join(pathlib.Path(__file__).parent.resolve(), 'data', f'{ref}.gbk')
    ref_features[ref] = {}
    if os.path.isfile(ref_file):
        for record in SeqIO.parse(ref_file, 'genbank'):
            k = '.'.join([
                os.path.splitext(os.path.basename(ref_file))[0],
                record.name,
            ])
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers and 'gene' not in feature.qualifiers:
                    feature.qualifiers['gene'] = deepcopy(feature.qualifiers['locus_tag'])
                if 'gene' in feature.qualifiers:
                    feature.hybranref = k
                    feature.references = {k: record.seq}
                    ref_features[ref][feature.qualifiers['gene'][0]] = feature


features = {
    'AZ20': {
        'secD': {
            'ratt': SeqFeature(
                FeatureLocation(3714209, 3716770, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['ECOLIN_02375'],'gene':['secD']}
            ),
        },
        'AZ20_03933': {
            'ratt': SeqFeature(
                FeatureLocation(3950866, 3951688, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['ECOLIN_01320'], 'gene':['ORF0033::ECOLIN_01320']},
            ),
            'prokka': SeqFeature(
                FeatureLocation(3950866, 3951688, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_03835'], 'gene':['ORF0033']},
            ),
        },
        'ECOLIN_18975': {
            'ratt': SeqFeature(
                FeatureLocation(333781, 336229, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['ECOLIN_18975'], 'gene':['ECOLIN_18975']},
            ),
        },
        'ECOLIN_18965': {
            'ratt': SeqFeature(
                FeatureLocation(334825, 336229, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['ECOLIN_18965'], 'gene':['ECOLIN_18965']},
            ),
        },
        'ECOLIN_25305': {
            'ratt': SeqFeature(
                FeatureLocation(622418, 623135, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['ECOLIN_25305'], 'gene':['ECOLIN_25305']},
            ),
        },
        'garD': {
            'ratt': SeqFeature(
                FeatureLocation(622418, 623990, strand=-1),
                type='CDS',
                qualifiers={
                    'locus_tag':['ECOLIN_17355'], 'gene':['garD'],
                    'note': [
                        "Hybran/Pseudoscan: Locus does not have reference-corresponding start | Poor blastp match at 95% identity and 95% coverage thresholds | No internal stop codons and ends with a valid stop codon | Locus divisible by three",
                        "D31 VS1 VE1 RCS0 RCE1 BOK0",
                    ],
                    'pseudo': [],
                },
            ),
        },
    },
    '1-0006': {
        'dnaA': {
            'ratt': SeqFeature(
                FeatureLocation(0, 1524, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0001'],'gene':['dnaA']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(33, 1524, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_00001'],'gene':['dnaA']}
            ),
        },
        'gyrB': {
            'abinit': SeqFeature(
                FeatureLocation(33, 1524, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_00007'],'gene_synonym':['gyrB1']}
            ),
        },
        'Rv0007': {
            'ratt': SeqFeature(
                FeatureLocation(11271, 12186, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0007'],'gene':['Rv0007']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(11880, 12186, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_00010'],'gene':['Rv0007'],'pseudo':[""]}
            ),
        },
        'galTb': {
	    'ratt': SeqFeature(
		FeatureLocation(712798, 713308, strand=1),
		type='CDS',
		qualifiers={'locus_tag':['Rv0619'],'gene':['galTb']}
	    ),
	},
        'Rv0061c': {
            'ratt': SeqFeature(
                FeatureLocation(66831, 67032, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0061c'],'gene':['Rv0061c']}
            ),
        },
        'Rv0074': {
            'ratt': SeqFeature(
                FeatureLocation(80866, 82276, strand=1),
                type='CDS',
                qualifiers={
                    'locus_tag':['Rv0074'],'gene':['Rv0074'],
                    'note': [
                        "Internal stop detected in the following codon(s): 119",
                    ],
                    'pseudo': [''],
                }
            ),
        },
        'Rv0071': {
            'ratt': SeqFeature(
                FeatureLocation(81235, 82276, strand=1),
                type='CDS',
                qualifiers={
                    'locus_tag':['Rv0071'],'gene':['Rv0071'],
                    'note': [
                        ("Hybran/Pseudoscan: Locus does not have reference-corresponding end | "
                         "Has a frameshift mutation leading to a delayed stop codon"),
                        "D31 VE1 RCS1 RCE0 BOK0",
                    ],
                    'pseudo': [''],
                }
            ),
        },
        'PPE5': {
            'ratt': SeqFeature(
                FeatureLocation(366754, 373353, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0304c'],'gene':['PPE5']}
            ),
        },
        'PPE6': {
            'ratt': SeqFeature(
                FeatureLocation(366753, 376299, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0305c'],'gene':['PPE6']}
            ),
        },
        'pip': {
            'ratt': SeqFeature(
                FeatureLocation(936376, 937360, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0840c'],'gene':['pip']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(936376, 937021, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_00896'],'gene':['pip']}
            ),
        },
        'Rv0907': {
            'ratt': SeqFeature(
                # Actually join(1012633..1012764,1012769..1014232)
                # query has a frameshift insertion, so combined interval isn't divisible by 3
                FeatureLocation(1012632, 1014232, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0907'],'gene':['Rv0907']}
                ),
            'abinit': SeqFeature(
                FeatureLocation(1012543, 1014232, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_00971'],'gene':['Rv0907']}
                ),
        },
        'Rv1075c': {
	    'ratt': SeqFeature(
		FeatureLocation(1202151, 1203042, strand=-1),
		type='CDS',
		qualifiers={'locus_tag':['Rv1075c'],'gene':['Rv1075c']}
	    ),
	},
        'PE10': {
	    'ratt': SeqFeature(
		FeatureLocation(1217413, 1217872, strand=1),
		type='CDS',
		qualifiers={'locus_tag':['Rv1089'],'gene':['PE10']}
	    ),
	},
        'Rv1150': {
            'ratt': SeqFeature(
                FeatureLocation(1280928, 1281480, strand=1),
                type='CDS',
                qualifiers={
                    'locus_tag':['L_01235'], 'gene':['Rv1150'], 'pseudo':[''],
                    'note': [
                        "Hybran/Pseudoscan: Reference gene is pseudo | Has reference-corresponding start and stop | Both this sequence and the reference's are divisible by three.",
                    ],
                }
            ),
            'prokka': SeqFeature(
                FeatureLocation(1280616, 1281480, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_01235'], 'gene':['Rv1041c']}
            ),
        },
        'PPE38': {
            'abinit': SeqFeature(
                FeatureLocation(2636045, 2637140, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_02521'],'gene':['PPE38']}
            ),
        },
        'PPE47': {
            'abinit': SeqFeature(
                FeatureLocation(3374234, 3375542, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_03201'],'gene':['PPE47']}
            ),
        },
        'Rv1453': {
            'ratt': SeqFeature(
                FeatureLocation(1642703, 1643969, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv1453'],'gene':['Rv1453']}
                ),
            'abinit': SeqFeature(
                FeatureLocation(1642670, 1643969, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_01557'],'gene':['Rv1453']}
                ),
        },
        'PE21': {
            'ratt': SeqFeature(
                FeatureLocation(2362218, 2362392, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv1548c'],'gene':['PE21'],'pseudo':['']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(2360919, 2362392, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_02249'],'gene':['PE21'],'pseudo':[''],
                            'note': [
                                "Has a frameshift mutation leading to a delayed stop codon",
                            ],
                }
            ),
        },
        'mce1R': {
	    'ratt': SeqFeature(
		FeatureLocation(192566, 193259, strand=-1),
		type='CDS',
		qualifiers={'locus_tag':['Rv0165c'],'gene':['mce1R']}
	    ),
	},
        'Rv1945': {
            'ratt': SeqFeature(
                FeatureLocation(2186568, 2187933, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv1945'],'gene':['Rv1945']}
                ),
            'abinit': SeqFeature(
                FeatureLocation(2186484, 2187933, strand=1),
                type='CDS',
                qualifiers={
                    'locus_tag':['L_02063'],'gene':['Rv1148c'],
                    'note': [
                        'Hybran/Pseudoscan: Strong blastp match at 95% identity and 95% coverage thresholds | Locus does not have reference-corresponding start',
                        'D31 VS1 VE1 RCS0 RCE1 BOK1',
                    ]
                }
                ),
        },
        'ORF0004': {
            'abinit': SeqFeature(
                # original was pseudo and FeatureLocation(2444022, 2444907, strand=-1)
                # before coord_check fixed it
                FeatureLocation(2444022, 2445009, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_02335'],'gene':['ORF0004']}
                ),
        },
        'Rv2079': {
	    'ratt': SeqFeature(
		FeatureLocation(2339465, 2341292, strand=1),
		type='CDS',
		qualifiers={'locus_tag':['Rv2079'],'gene':['Rv2079']}
	    ),
	},
        'Rv2081c': {
	    'ratt': SeqFeature(
		FeatureLocation(2342287, 2342617, strand=-1),
		type='CDS',
		qualifiers={'locus_tag':['Rv2081c'],'gene':['Rv2081c']}
	    ),
	},
        'Rv2134c': {
	    'ratt': SeqFeature(
		FeatureLocation(2397623, 2398106, strand=-1),
		type='CDS',
		qualifiers={'locus_tag':['Rv2134c'],'gene':['Rv2134c']}
	    ),
	},
        'Rv2180c': {
            'ratt': SeqFeature(
                FeatureLocation(2444022, 2445033, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv2180c'],'gene':['Rv2180c']}
                ),
        },
        'Rv1877': {
            'ratt': SeqFeature(
                FeatureLocation(2112334, 2114834, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv1877'],'gene':['Rv1877']}
            )
        },
        'Rv2280': {
	    'ratt': SeqFeature(
		FeatureLocation(2553287, 2553950, strand=1),
		type='CDS',
		qualifiers={'locus_tag':['Rv2280'],'gene':['Rv2280']}
	    ),
	},
        'Rv2437': {
	    'ratt': SeqFeature(
		FeatureLocation(2735341, 2736310, strand=1),
		type='CDS',
		qualifiers={'locus_tag':['Rv2437'],'gene':['Rv2437']}
	    ),
	},
        'Rv2879c': {
            'ratt': SeqFeature(
                FeatureLocation(3182302, 3183397, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv2879c'],'gene':['Rv2879c']}
            )
        },
        'Rv2880c': {
            'ratt': SeqFeature(
                FeatureLocation(3182302, 3183397, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv2880c'],'gene':['Rv2880c']}
            )
        },
        'Rv3181c': {
            'abinit': SeqFeature(
                FeatureLocation(3548089, 3548338, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_03370'],'gene':['Rv3181c'],
                            'pseudo':['']}
            )
        },
        'Rv3327': {
            'ratt': SeqFeature(
                FeatureLocation(3707086, 3709176, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['Rv3327'], 'gene':['Rv3327']}
            )
        },
        'PE_PGRS50': {
            'ratt_raw': SeqFeature(
                FeatureLocation(3741108, 3746955, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv3345c'],'gene':['PE_PGRS50']}
            ),
            'final': SeqFeature(
                FeatureLocation(3741108, 3746955, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv3345c'],'gene':['PE_PGRS50::PE_PGRS49']}
            ),
        },
        'Rv3777': {
            'abinit': SeqFeature(
                FeatureLocation(4230776, 4231754, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_04009'],'gene':['Rv3777']}
            )
        },
        'Rv3785': {
            'abinit': SeqFeature(
                FeatureLocation(4239976, 4240378, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_04018'],'gene':['Rv3785']}
            )
        },
    },
    '1-0009': {
        'PPE54': {
            'ratt': SeqFeature(
                FeatureLocation(3730264, 3741334, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['PPE54'],'gene':['PPE54']}
            )
        },
    },
    '1-0047': {
        'PE_PGRS54': {
            'abinit': SeqFeature(
                FeatureLocation(3959135, 3959759, strand=1),
                type='CDS',
                qualifiers={
                    'locus_tag':['L_03743'],
                    'translation':['MASEGGAGGQGGDGGDGGEGGGAGFGSGVAGAAGAGGNGGKGGDGGTGGTGGTNFAGGQGGAGGRGGAGGNGANGVGDNAAGGDGGNGGAGGLGGGGGTGGTNGNGGLGGGGGNGGAGGAGGTPTGSGTEGTGGDGGDAGAGGNGGSATGVGNGGNGGDGGNGGDGGNGAPGGFGGGAGAGGLGGSGAGGGTDGDDGNGGSPGTDGS'],
                }
            ),
        },
    },
    '2-0031': {
        'PPE34': {
            'abinit': SeqFeature(
                FeatureLocation(1927378, 1927894, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_01804'],'gene':['PPE34']}
            )
        },
    },
    '4-0041': {
        'rrf': {
            'ratt': SeqFeature(
                FeatureLocation(2935204, 2935319, strand=-1),
                type='rRNA',
                qualifiers={'locus_tag':['MTB000021'],'gene':['rrf']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(2935204, 2935316, strand=-1),
                type='rRNA',
                qualifiers={'locus_tag':['L_02791']}
            ),
        },
        'Rv1718': {
            'ratt': SeqFeature(
                FeatureLocation(2476948, 2477773, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv1718'],'gene':['Rv1718']}
            ),
            'abinit': SeqFeature(
                FeatureLocation(2476948, 2477767, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['L_02383'],'gene':['Rv1718']}
            ),
        },
    },
}

abinit_blast_results = {
    '1-0006': {
        'L_00001': {'iden': 99.8, 'qcov': 100.0, 'scov': 97.8},
        'L_00010': {'iden': 99.01, 'qcov': 100.0, 'scov': 33.224}, # Rv0007
        'L_01557': {'iden': 99.525, 'qcov': 97.454, 'scov': 100.0}, # Rv1453
        'L_02063': {'iden': 96.888, 'qcov': 100.0, 'scov': 100.0}, # Rv1148c
        'L_02335': {'iden': 100.0, 'qcov': 100.0, 'scov': 100.0}, # ORF0004 (Rv3475)
    },
    '4-0041': {
        'L_02383': {'iden': 100.0, 'qcov': 100.0, 'scov': 100.0},
    },
}

ratt_blast_results = {
    '1-0006': {
        'Rv0001': {'iden': 99.8, 'qcov': 100.0, 'scov': 100.0},
        'Rv0007': {'iden': 99.671, 'qcov': 100.0, 'scov': 100.0},
        'Rv1453': {'iden': 99.525, 'qcov': 100.0, 'scov': 100.0},
        'Rv2180c': {'iden': 41.176, 'qcov': 4.762, 'scov': 5.424}, # overlaps with L_02235
        'Rv2946c': {'iden': 99.876, 'qcov': 100.0, 'scov': 100.0},

    },
    '4-0041': {
        'Rv1718': {'iden': 100.0, 'qcov': 99.270, 'scov': 100.0},
    },
}
