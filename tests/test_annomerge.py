from copy import deepcopy

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
                                             'pseudo' in feature_list[2].qualifiers.keys() and feature_list[2].qualifiers['pseudo'] == [''] # Rv0739


@pytest.mark.parametrize('feature', [
    'Rv0205_1',
    'Rv0205_2',
    'rplB',
    'mamB',
])
def test_validate_prokka_feature_annotation(feature):
    inputs = dict(
        # Rv0205 in this isolate is disrupted.
        # These two gene fragments should be tagged with /pseudo qualifiers by this function and
        # later collected under a single /gene entry.
        Rv0205_1 = SeqFeature(
            FeatureLocation(ExactPosition(245198), ExactPosition(245834), strand=1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_00229"], gene=["Rv0205"], codon_start=['1'], transl_table='11',
                translation=['MSASLDDASVAPLVRKTAAWAWRFLVILAAMVALLWVLNKFEVIVVPVLLALMLSALLVPPVDWLDSRGLPHAVAVTLVLLSGFAVLGGILTFVVSQFIAGLPHLVTEVERSIDSARRWLIEGPAHLRGEQIDNAGNAAIEALRNNQAKLTSGALSTAATITELVTAAVLILFTLIFFLYGAGASGSTSRRPSRPASVTECVRRGAPVMRR']
            )
        ),
        Rv0205_2 = SeqFeature(
            FeatureLocation(ExactPosition(245797), ExactPosition(246301), strand=1), type='CDS',
            qualifiers = dict(
                locus_tag=["L_00230"], gene=["Rv0205"], codon_start=['1'], transl_table='11',
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

    expected_feature = deepcopy(inputs[feature])
    if feature in ['Rv0205_1', 'Rv0205_2']:
        expected_feature.qualifiers['pseudo'] = ['']

    expected = dict(
        Rv0205_1 = (expected_feature, False, ''),
        Rv0205_2 = (expected_feature, False, ''),
        rplB = (expected_feature, False, ''),
        mamB = (expected_feature, False, ''),
    )

    prokka_noref = dict()
    ratt_blast_results = dict()
    reference_locus_list = []
    seq_ident = seq_covg = 95
    ref_temp_fasta_dict = dict(
        Rv0205='gene-seqs/Rv0205.fasta',
        Rv0704='gene-seqs/Rv0704.fasta',
        Rv2024c='gene-seqs/Rv2024c.fasta',
    )
    reference_gene_locus_dict = dict(
        Rv0205='Rv0205',
        rplB='Rv0704',
        mamB='Rv2024c',
    )
    reference_locus_gene_dict = dict(
        Rv0205='Rv0205',
        Rv0704='rplB',
        Rv2024c='mamB',
    )

    annomerge.prokka_blast_list = []
    results = annomerge.validate_prokka_feature_annotation(
        deepcopy(inputs[feature]),
        prokka_noref,
        reference_gene_locus_dict,
        reference_locus_gene_dict,
        ref_temp_fasta_dict,
        ratt_blast_results,
        reference_locus_list,
        seq_ident,
        seq_covg
    )
    # __eq__ is not implemented for the SeqFeature object itself,
    # so check the qualifiers and then check the rest of the return values
    assert results[0].qualifiers == expected_feature.qualifiers \
        and results[1:] == expected[feature][1:]
