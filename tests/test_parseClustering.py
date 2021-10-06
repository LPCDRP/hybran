from io import StringIO

from hybran import parseClustering


# sample input for prepare_for_eggnog

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
                                       target_genomes = ['a'],
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
