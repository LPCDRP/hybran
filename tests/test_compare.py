from Bio.SeqFeature import SeqFeature
import pytest

from hybran import compare


@pytest.mark.parametrize('case', [
    'two_word_code',
    'multiple_codes',
    'not_relevant',
])
def test_pgap_np(case):
    note = {
        'two_word_code': (
            "internal stop; Derived by automated computational "
            "analysis using gene prediction method: Protein Homology. "
            "GO_function: GO:0003677 - DNA binding [Evidence IEA]; "
            "GO_function: GO:0004519 - endonuclease activity [Evidence "
            "IEA]"
        ),
        'multiple_codes': (
            "internal stop; incomplete; partial on complete "
            "genome; missing C-terminus; Derived by automated "
            "computational analysis using gene prediction method: "
            "Protein Homology."
        ),
        'not_relevant': (
            "Essential for efficient processing of 16S rRNA; "
            "Derived by automated computational analysis using gene "
            "prediction method: Protein Homology. "
            "GO_function: GO:0003723 - RNA binding [Evidence IEA]; "
            "GO_process: GO:0006364 - rRNA processing [Evidence IEA]"
        ),
    }
    expected = {
        'two_word_code': 'internal_stop',
        'multiple_codes': 'internal_stop;incomplete',
        'not_relevant': '.',
    }

    assert compare.pgap_np(
        SeqFeature(qualifiers={'note': [ note[case] ]})
    ) == expected[case]
