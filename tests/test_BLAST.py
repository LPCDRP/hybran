from Bio import SeqIO

from hybran import BLAST


def test_blastp():
    record = SeqIO.read("gene-seqs/Rv2434c.fasta", 'fasta')
    assert BLAST.blastp(record, "gene-seqs/Rv2434c.fasta", 95, 95) == \
        ['Rv2434c\tRv2434c\t100.000\t481\t0\t0\t1\t481\t1\t481\t0.0\t966\t481\t481\t100.0\t100.0']
