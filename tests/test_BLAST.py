from collections import defaultdict
import os

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord
import pytest

from hybran import BLAST
from hybran.bio import SeqIO

from .data_features import *


@pytest.mark.skipif(not os.path.isfile("data/H37Rv.faa"), reason="reference proteome data file not available")
@pytest.mark.parametrize('feature', [
    'Rv0205_1',
    'rplB',
    'mamB',
    'no_hit',
])
def test_reference_match(feature):
    inputs = dict(
        # Rv0205 in this isolate is disrupted.
        Rv0205_1 = SeqFeature(
            FeatureLocation(ExactPosition(245198), ExactPosition(245834), strand=1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_00229"], gene=["Rv0205"], codon_start=['1'], transl_table='11',
                pseudo=[''],
                translation=['MSASLDDASVAPLVRKTAAWAWRFLVILAAMVALLWVLNKFEVIVVPVLLALMLSALLVPPVDWLDSRGLPHAVAVTLVLLSGFAVLGGILTFVVSQFIAGLPHLVTEVERSIDSARRWLIEGPAHLRGEQIDNAGNAAIEALRNNQAKLTSGALSTAATITELVTAAVLILFTLIFFLYGAGASGSTSRRPSRPASVTECVRRGAPVMRR']
            )
        ),
        Rv0205_2 = SeqFeature(
            FeatureLocation(ExactPosition(245797), ExactPosition(246301), strand=1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_00230"], gene=["Rv0205"], codon_start=['1'], transl_table='11',
                pseudo=[''],
                translation=['MRAAGRAGYASLIGYARATFLVALTDAAGVGAGLAVMGVPLALPLASLVFFGAFIPLIGAVVAGFLAVVVALLAKGIGYALITVGLLIAVNQLEAHLLQPLVMGRAVSIHPLAVVLAIAAGGVLAGVVGALLAVPTVAFFNNAVQVLLGGNPFADVADVSSDHLTEV']
            )
        ),
        # This is a normal gene
        rplB = SeqFeature(
            FeatureLocation(ExactPosition(809202), ExactPosition(810045), strand=1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_00759"], gene=["rplB"], codon_start=['1'], transl_table='11',
                translation=['MAIRKYKPTTPGRRGASVSDFAEITRSTPEKSLVRPLHGRGGRNAHGRITTRHKGGGHKRAYRMIDFRRNDKDGVNAKVAHIEYDPNRTARIALLHYLDGEKRYIIAPNGLSQGDVVESGANADIKPGNNLPLRNIPAGTLIHAVELRPGGGAKLARSAGSSIQLLGKEASYASLRMPSGEIRRVDVRCRATVGEVGNAEQANINWGKAGRMRWKGKRPSVRGVVMNPVDHPHGGGEGKTSGGRHPVSPWGKPEGRTRNANKSSNKFIVRRRRTGKKHSR']
            )
        ),
        # This gene is truncated in the reference.
        mamB = SeqFeature(
            FeatureLocation(ExactPosition(2272132), ExactPosition(2276953), strand=-1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_02174"], gene=["mamB"], codon_start=['1'], transl_table='11',
                translation=['MGSVHDVIEAFRKAPSNAERGTKFEQLMVRYFELDPTMAQQYDAVWRWIDWPERRGRTDTGIDLVARERDTGNYTAIQCKFYEPTHTLAKGDIDSFFTASGKTGFTNRVIISTTDRWGRNAEDALADQLVPVQRIGMAEIAESPIDWDIAWPAGDLQVNLTPAKRHELRPHQQQAIDAVFRGFAVGNDRGKLIMACGTGKTFTALKIAERIAADNGGSARILLLVPSISLLSQTLREWTAQSELDVRAFAVCSDTKVSRSAEDYHVHDVPIPVTTDARVLLHEMAHRRRAQGLTVVFCTYQSLPTVAKAQRLGVDEFDLVMCDEAHRTTGVTLAGDDESNFVRVHDGQYLKAARRLYMTATPRIFTESIKDRADQHSAELVSMDDELTFGPEFHRLSFGEAVERGLLTDYKVMVLTVDQGVIAPRLQQELSGVSGELMLDDASKIVGCWNGLAKRSGTGIVAGEPPMRRAVAFAKDIKTSKQVAELFPKVVEAYRELVDDGPGLACSVRHVDGTFNALVRNEQLAWLKGVVAEDECRILSNARCLSEGVDVPALDAVLFLNPRNSIVDVVQSVGRVMRKSPGKDYGYVILPVAVPEGVEPSAALADNKRFKVVWQVLNALRSHDERFDAMVNSIALNVKPTKTGEGSDKLLGGHIGPTSDEAGPAVAEQLAMFSLSQWQEAIYARIVDKVGTRTYWEQWAADVADIAATLTTRIHALLGGADATAAAAFEQFLAGLRDNLNDSITPDDAISMLSQHLITKPVFDALFAGHDFASHNPVSRAMQKMVDTVGGAGLEAETARLEGFYESVRRRAGEVTSAEGKQQVIAELYEKFFRIGFKKQAEALGIVYTPVEVVDFIVRAADFVSRKHFGRGLTDEGVHILDGFAGTGTFITRLLQSDLITAADLTRKYSQELHANEIMLLAYYIAAVNIESTYHALAGKTADADAYEPFPGMALADTFQISEAGDSMDAIMFPYNNARILRQLATPISVIIGNPPYSVGQSSANDLNANVKYPTLDGRIEQTYAKRSTAQLKNSLYDSYIRAFRWATDRIGDNGVVGFVSNGGYIDGNTADGMRLSLADDYAAVYVYNLRGNQRTAGELSRQEGGKVFGGGSRNTVAIFLGIKDPKHSGPCDVLYRDIGDYLSREEKLRIVGDGYLDTVEWQTVTPNLHGDWVNQRDDAFSAWPVIGDKKAALDVTRVFANYSAGLKTSRDAWCYNFSRGALEANIGRTIDFYNSEVDRINEIRGRDAKTPPVDALITVDSAKFSWDRINKRQVAQGIRIEFAPAGMRLGTYRPFTKEHAYLDPNQQLNNCTYQLPSMFPTPEHGNVGYYVVGMGSDKPFSCLMLNAIPDLAFWGSSNGQFFPRWTYEKTEPRDGELDFESTTNAEVDDHGYRRVDNITGVILKLYRDTIGDQVTKDDIFYYVYGLLHDPAYRTKYAADLKKMLPHIPTPETRERFDQLASAGRKLADLHVGYESVKPYPLDVQLKPGADPEDRETWRVEKMKWKSKQDHSTIIYNSRVTIAGIPDEAERYLLGSRSALGWIIDRYRVTTDKASGIVNDPNDWCDEHANPTYIVDLIKKVTTVSVETMKIVDSIVALASAGSDST']
            )
        ),
        no_hit = features['1-0047']['PE_PGRS54']['abinit'],
    )

    expected = {
        'Rv0205_1': (
            None, False,
            {
                'Rv0205': {
                    'iden': 92.611,
                    'qcov': 96.2085308056872,
                    'scov': 55.04087193460491,
                },
                'Rv2018': {
                    'iden': 32.143,
                    'qcov': 26.540284360189574,
                    'scov': 20.502092050209207,
                },
                'Rv3736': {
                    'iden': 40.0,
                    'qcov': 14.218009478672986,
                    'scov': 9.91501416430595,
                },
                'asnB': {
                    'iden': 22.656,
                    'qcov': 60.66350710900474,
                    'scov': 17.177914110429448,
                },
                'betP': {
                    'iden': 28.814,
                    'qcov': 24.644549763033176,
                    'scov': 9.949409780775717,
                },
                'bkdC': {
                    'iden': 43.75,
                    'qcov': 15.165876777251185,
                    'scov': 6.106870229007633,
                },
                'gnd1': {
                    'iden': 47.619,
                    'qcov': 9.95260663507109,
                    'scov': 4.329896907216495,
                },
                'pks3': {
                    'iden': 29.73,
                    'qcov': 34.12322274881517,
                    'scov': 15.163934426229508,
                },
            },
        ),
        'rplB': (
            'rplB', False,
            {'rplB': {'iden': 100.0, 'qcov': 100.0, 'scov': 100.0}},
        ),
        'mamB': (
            'mamB', True,
            {
                'mamB': {
                    'iden': 99.605,
                    'qcov': 31.506849315068493,
                    'scov': 98.25242718446601,
                },
            },
        ),
        'no_hit': (
            None, False,
            defaultdict(lambda : {
                'iden':0.,
                'qcov':0.,
                'scov':0.,
            })
        ),
    }

    for ex_hit, ex_lowcovg, ex_results in expected.values():
        for key in ex_results:
            for datum in ['iden', 'qcov', 'scov']:
                if ex_results[key][datum] not in [0., 100.]:
                    ex_results[key][datum] = pytest.approx(ex_results[key][datum])

    (result, low_covg, hits) = BLAST.reference_match(
        query=SeqRecord(Seq(inputs[feature].qualifiers['translation'][0])),
        subject="data/H37Rv.faa",
        seq_ident=95,
        seq_covg=95,
    )
    hits = dict(hits)
    assert (result, low_covg, hits) == expected[feature]

def test_summarize():
    blast_results = \
        ['Rv2434c\tRv2434c\t100.000\t481\t0\t0\t1\t481\t1\t481\t0.0\t966\t481\t481\t100.0\t100.0',
         'Rv0000\tRv0000\t100.000\t481\t0\t0\t1\t481\t1\t481\t0.0\t966\t481\t481\t100.0\t100.0',
         'Rv2434c\tRv2434c\t96.000\t481\t0\t0\t1\t481\t1\t481\t0.0\t966\t481\t481\t100.0\t100.0',
         ]
    assert BLAST.summarize(blast_results) == \
        {'Rv2434c': {
            'iden': 100.0,
            'scov': 100.0,
            'qcov': 100.0,
         },
         'Rv0000': {
             'iden': 100.0,
             'scov': 100.0,
             'qcov': 100.0,
         }
        }

def test_blastp():
    record = SeqIO.read("gene-seqs/Rv2434c.fasta", 'fasta')
    assert BLAST.blastp(record, "gene-seqs/Rv2434c.fasta", 95, 95) == \
        (['Rv2434c\tRv2434c\t100.000\t481\t0\t0\t1\t481\t1\t481\t0.0\t966\t481\t481\t100.0\t100.0'],
         ['Rv2434c\tRv2434c\t0\t0\t0\t0\t\t\t\t\t1\t0\t1\t1\t0.0\t0.0'],
         [])
