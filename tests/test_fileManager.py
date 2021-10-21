import os
import tempfile

import pytest

from hybran import fileManager


@pytest.fixture(scope="session")
def refdir(tmp_path_factory):
    refd = tmp_path_factory.mktemp("refs")
    fofn = refd.joinpath('refs.fofn')
    fullpaths = []
    for f in ['ref1.fasta','ref2.fna', 'ref3.fa','garbage.txt']:
        reffile = refd.joinpath(f)
        reffile.touch()
        if f != 'garbage.txt':
            fullpaths.append(str(reffile.resolve()))
    fofn.write_text('\n'.join(fullpaths))
    return refd


@pytest.mark.parametrize("test_input", [
    ['.'],
    ['ref1.fasta', 'ref2.fna', 'ref3.fa'],
    ['refs.fofn'],
])
def test_genome_list(test_input, refdir):
    os.chdir(refdir)
    assert sorted(fileManager.genome_list(test_input)) == \
        [str(refdir.joinpath(_)) for _ in ['ref1.fasta','ref2.fna','ref3.fa']]
