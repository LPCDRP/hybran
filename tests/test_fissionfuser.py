from collections import defaultdict
import os

from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    ExactPosition,
)
import pytest

from hybran import (
    annomerge,
    config,
    fissionfuser,
)
from hybran.bio import (
    SeqIO,
    translate,
)
from hybran.util import keydefaultdict

from .data_features import *


@pytest.mark.parametrize('gene_list', [
    'complementary_fragments',
    'complementary_fragments_one_unnamed',
    'fails_final_coord_check_inframe_overlap',
#    'seemingly_complete_fragment',
    'independent_copies',
])
@pytest.mark.skipif(not os.path.isfile("data/1-0006.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/2-0031.fasta"), reason="test genome sequence not available")
@pytest.mark.skipif(not os.path.isfile("data/H37Rv.fasta"), reason="test reference genome sequence not available")
def test_fissionfuser(gene_list, tmp_path):
    inputs = {
        # post coordinate correction, so they both have the same start position
        'complementary_fragments' : [
            SeqFeature(
                FeatureLocation(ExactPosition(2275540), ExactPosition(2277261), strand=-1), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02173'],
                    'pseudo': [''],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(2276448), ExactPosition(2277261), strand=-1), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02174'],
                    'pseudo': [''],
                }
            ),
        ],
        'complementary_fragments_one_unnamed': [
            SeqFeature(
                FeatureLocation(ExactPosition(1104310), ExactPosition(1104607), strand=1), type='CDS', qualifiers={
                    'locus_tag': ['L_01053'],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(1104558), ExactPosition(1105056), strand=1), type='CDS', qualifiers={
                    'gene': ['Rv0986'],
                    'locus_tag': ['L_01054'],
                    'pseudo': [''],
                }
            ),
        ],
        'fails_final_coord_check_inframe_overlap' : [
            SeqFeature(
                # pre-correction: FeatureLocation(ExactPosition(3522323), ExactPosition(3523235), strand=-1)
                FeatureLocation(ExactPosition(3522323), ExactPosition(3523418), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03351'],
                    'pseudo': [''],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(3523256), ExactPosition(3523418), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03352'],
                    'pseudo': [''],
                }
            ),
        ],
        'seemingly_complete_fragment': [
        ],
        'independent_copies': [
            SeqFeature(
                FeatureLocation(ExactPosition(1483334), ExactPosition(1484933), strand=-1), type='CDS', qualifiers={
                    'locus_tag': ['L_01415'],
                }
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(1484939), ExactPosition(1486616), strand=-1), type='CDS', qualifiers={
                    'gene': ['Rv1319c'],
                    'locus_tag': ['L_01416'],
                }
            ),
        ],
    }
    source_genomes = {
        'complementary_fragments':'1-0006',
        'complementary_fragments_one_unnamed':'1-0006',
        'fails_final_coord_check_inframe_overlap':'2-0031',
        'independent_copies':'1-0006',
    }
    config.hybran_tmp_dir = tmp_path
    ref_genome = defaultdict(lambda :'H37Rv')
    source_genome = source_genomes[gene_list]
    record_sequence = list(SeqIO.parse(f'data/{source_genome}.fasta', 'fasta'))[0]
    config.cnf.genetic_code = 11
    ref_annotation = keydefaultdict(annomerge.ref_fuse)
    ref_annotation.update(ref_features[ref_genome[gene_list]])

    for f in inputs[gene_list]:
        f.source = 'H37Rv.NC_000962.3' # TODO - clean this up
        f.location.ref = record_sequence.id
        f.references = {record_sequence.id: record_sequence.seq}
        f.qualifiers['translation'] = [ str(translate(
            f.extract(record_sequence.seq),
            table=config.cnf.genetic_code,
            to_stop=True,
        )) ]


    expected = {
        'complementary_fragments': ([
            SeqFeature(
                FeatureLocation(ExactPosition(2275540), ExactPosition(2277261), strand=-1, ref='1'), type='CDS', qualifiers={
                    'gene': ['dosT'],
                    'locus_tag': ['L_02174'],
                    'note': [
                        ('Hybran/Pseudoscan:description:Internal stop detected at '
                         'codon(s) 270 283 353 380 | Locus has invalid reading frame-- not '
                         'divisible by three | Locus has reference-corresponding start and '
                         'end | Poor blastp match at 50 bitscore and 80% coverage thresholds '
                         '| Locus is 1 base pair(s) shorter than the reference'),
                        'Hybran/Pseudoscan:evidence:not_div_by_3;internal_stop',
                        'Hybran/Pseudoscan:barcode:D30;VS1;VE0;RCS1;RCE1;BOK0',
                    ],
                    'pseudo': [''],
                    'translation': ['MTHPDRANVNPGSPPLRETLSQLRLRELLLEVQDRIEQIVEGRDRLDGLIDAILAITSGLKLDATLRAIVHTAAELVDARYGALGVRGYDHRLVEFVYEGIDEETRHLIGSLPEGRGVLGALIEEPKPIRLDDISRHPASVGFPLHHPPMRTFLGVPVRIRDEVFGNLYLTEKADGQPFSDDDEVLVQALAAAAGIAVDNARLFEESRTREAWIEATRDIGTQMLAGADPAMVFRLIAEEALTLMAGAATLVAVPLDDKRRLARSTTWSS'],
                }
            )],
            [{
                'feature':inputs['complementary_fragments'][0],
                'evid':'overlapping_inframe',
                'remark':'combined fission fragments',
            }],
        ),
        'complementary_fragments_one_unnamed': ([
            SeqFeature(
                FeatureLocation(ExactPosition(1104310), ExactPosition(1105056), strand=1, ref='1'), type='CDS', qualifiers={
                    'gene': ['Rv0986'],
                    'locus_tag': ['L_01054'],
                    'note': [
                        ('Hybran/Pseudoscan:description:Internal stop detected at '
                         'codon(s) 98 168 199 208 215 223 227 | Locus has invalid reading '
                         'frame-- not divisible by three | Locus has reference-corresponding '
                         'start and end | Poor blastp match at 50 bitscore and 80% coverage '
                         'thresholds | Locus is 1 base pair(s) shorter than the reference'),
                        'Hybran/Pseudoscan:evidence:not_div_by_3;internal_stop',
                        'Hybran/Pseudoscan:barcode:D30;VS1;VE0;RCS1;RCE1;BOK0',
                    ],
                    'pseudo': [''],
                    'translation': ['MFRRDQIGIVFQFFNLIPTLTVLENITLPQELAGVSQRKAAVVARDLLEKVGMADRERTFPDKLSGGEQQRVAISRALAHNPMLVLADEPTGNLDSDTGDKVLDVLLDLTRQAGKTLIMATHSPSMTQHADRVVNLQGGRLIPAVNRENQTDQPASTILLPTSYE'],
                }
            )],
            [{
                'feature':inputs['complementary_fragments_one_unnamed'][0],
                'evid':'complementary_fragments',
                'remark':'combined fission fragments',
            }],
        ),
        'fails_final_coord_check_inframe_overlap': ([
            SeqFeature(
                FeatureLocation(ExactPosition(3522323), ExactPosition(3523418), strand=-1, ref='1'), type='CDS', qualifiers={
                    'gene': ['Rv3327'],
                    'locus_tag': ['L_03352'],
                    'note': [
                        ('Hybran/Pseudoscan:description:Internal stop detected at '
                         'codon(s) 53 | Locus has valid reading frame | Locus does not have '
                         'reference-corresponding end | Poor blastp match at 50 bitscore '
                         'and 80% coverage thresholds'),
                        'Hybran/Pseudoscan:evidence:internal_stop',
                        'Hybran/Pseudoscan:barcode:D31;VS1;VE0;RCS1;RCE0;BOK0',
                    ],
                    'pseudo': [''],
                    'translation': ['MVVVGTDAHKYSHTFVATDEVGRQLGEKTVKATTAGHATAIMWAREQFGLELI'],
                }
            )],
            [{
                'feature':inputs['fails_final_coord_check_inframe_overlap'][0],
                'evid':'overlapping_inframe',
                'remark':"combined fission fragments",
            }],
        ),
        'independent_copies': (
            inputs['independent_copies'],
            [],
        ),
    }

    # set the rival feature as the passing one. for these test cases, we're always looking at pairs
    if expected[gene_list][1]:
        expected[gene_list][1][0]['superior'] = expected[gene_list][0][0]

    assert fissionfuser.fissionfuser(
        inputs[gene_list],
        ref_annotation,
    ) == expected[gene_list]
