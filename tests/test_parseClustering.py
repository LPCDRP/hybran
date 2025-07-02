from io import StringIO

from hybran import parseClustering


def test_single_gene_clusters():
    single_gene_dict = {
        'H37Rv,Rv1076,LipU': [
            ['H37Rv', 'Rv1076', 'LipU'],
            ['1-0006', 'Rv1076', 'LipU']
        ],
        'H37Rv,Rv3763,LpqH': [
            ['H37Rv', 'Rv3763', 'LpqH'],
            ['1-0006', 'Rv3763', 'LpqH']
        ],
        'H37Rv,Rv2109c,Rv2109c': [
            ['H37Rv', 'Rv2109c', 'Rv2109c'],
            ['1-0006', '1-0006_02258', 'L_02258']
        ]
    }
    # initialize parseClustering's global variable
    parseClustering.isolate_update_dictionary = dict()

    parseClustering.single_gene_clusters(single_gene_dict)

    print(parseClustering.isolate_update_dictionary)
    assert parseClustering.isolate_update_dictionary == \
        {'1-0006':{'1-0006_02258': {'name':'Rv2109c', 'pseudo':False}}}


def test_prepare_for_eggnog():
    unannotated_seqs = StringIO(
""">MTB0084 False
MFLGQATIWAGWAILCGRLPVATGLAVFVGIWFVASP
>MTB0085 False
MTTLPALGYLSGALLHALVDYPTFSDRVLRCGAC
>MTB0086 True
MSLCRFGFQLDQSNLKLDQSNLTCKRISSIFTMV
>MTB0087 False
MVEMLGVCVVATDLDGHRIVLRDSEKTLKMRNR
>MTB0088 True
MRVLVTKPDGTQVEVHLDQGFRFLGTETVDND
"""
    )
    outfile = StringIO()
    parseClustering.prepare_for_eggnog(unannotated_seqs = unannotated_seqs,
                                       outfile = outfile)
    assert StringIO.getvalue(outfile) == \
""">MTB0084 False
MFLGQATIWAGWAILCGRLPVATGLAVFVGIWFVASP
>MTB0085 False
MTTLPALGYLSGALLHALVDYPTFSDRVLRCGAC
>MTB0087 False
MVEMLGVCVVATDLDGHRIVLRDSEKTLKMRNR
"""

if __name__=='__main__':
    test_prepare_for_eggnog()
