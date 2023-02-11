from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, CompoundLocation


ref_features = {
    'H37Rv': {
        'dnaA': SeqFeature(
            FeatureLocation(0, 1524, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv0001'],'gene':['dnaA']}
        ),
        'rrf': SeqFeature(
            FeatureLocation(1476898, 1477013, strand=1),
            type='rRNA',
            qualifiers={'locus_tag':['MTB000021'],'gene':['rrf']}
        ),
        'Rv1453': SeqFeature(
            FeatureLocation(1638380, 1639646, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv1453'],'gene':['Rv1453']}
        ),
        'Rv1718': SeqFeature(
            FeatureLocation(1944808, 1945627, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv1718'],'gene':['Rv1718']}
        ),
        'Rv3181c': SeqFeature(
            FeatureLocation(3549690,3550143 , strand=-1),
            type='CDS',
            qualifiers={'locus_tag':['Rv3181c'],'gene':['Rv3181c']}
        ),
        'Rv3777': SeqFeature(
            FeatureLocation(4222693, 4223680, strand=1),
            type='CDS',
            qualifiers={'locus_tag':['Rv3777'],'gene':['Rv3777']}
        ),
    },
}

features = {
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
                qualifiers={'locus_tag':['L_00007'],'gene':['gyrB1']}
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
                FeatureLocation(366754, 376299, strand=-1),
                type='CDS',
                qualifiers={'locus_tag':['Rv0305c'],'gene':['PPE6']}
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
        'Rv3777': {
            'abinit': SeqFeature(
                FeatureLocation(4230776, 4231754, strand=1),
                type='CDS',
                qualifiers={'locus_tag':['L_04009'],'gene':['Rv3777']}
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
        'L_01557': {'iden': 99.525, 'qcov': 97.454, 'scov': 100.0}, # Rv1453
    },
    '4-0041': {
        'L_02383': {'iden': 100.0, 'qcov': 100.0, 'scov': 100.0},
    },
}

ratt_blast_results = {
    '1-0006': {
        'Rv0001': {'iden': 99.8, 'qcov': 100.0, 'scov': 100.0},
        'Rv1453': {'iden': 99.525, 'qcov': 100.0, 'scov': 100.0},
    },
    '4-0041': {
        'Rv1718': {'iden': 100.0, 'qcov': 99.270, 'scov': 100.0},
    },
}
