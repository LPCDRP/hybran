from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, CompoundLocation


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
    },
}
