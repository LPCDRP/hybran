import os

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord
import pytest

from hybran import BLAST


@pytest.mark.skipif(not os.path.isfile("data/H37Rv.faa"), reason="reference proteome data file not available")
@pytest.mark.parametrize('feature', [
    'Rv0205_1',
    'rplB',
    'mamB',
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
        )
    )

    expected = {
        'Rv0205_1': (
            None, False,
            {
                'Rv0205': {'iden': 92.611, 'qcov': 60.23738872403561, 'scov': 55.04087193460491},
                'betP': {'iden': 28.814, 'qcov': 15.43026706231454, 'scov': 9.949409780775717},
                'pks3': {'iden': 32.308, 'qcov': 18.694362017804153, 'scov': 13.319672131147541},
            },
        ),
        'rplB': (
            'rplB', False,
            {'rplB': {'iden': 100.0, 'qcov': 71.06598984771574, 'scov': 100.0}},
        ),
        'mamB': (
            'mamB', False,
            {'mamB': {'iden': 99.605, 'qcov': 29.418604651162788, 'scov': 98.25242718446601}},
        ),
    }

    (result, pseudo, hits) = BLAST.reference_match(
        query=SeqRecord(inputs[feature]),
        subject="data/H37Rv.faa",
        seq_ident=95,
        seq_covg=95,
    )
    hits = dict(hits)
    assert (result, pseudo, hits) == expected[feature]

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
