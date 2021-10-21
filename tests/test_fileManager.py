import os
import tempfile

import pytest

from hybran import fileManager


@pytest.fixture(scope="session")
def refdir(tmp_path_factory):
    refd = tmp_path_factory.mktemp("refs")
    seq_fofn = refd.joinpath('refs.fofn')
    gbk_fofn = refd.joinpath('anns.fofn')
    fastapaths = []
    annpaths = []
    for f in ['ref1.fasta','ref2.fna', 'ref3.fa',
              'ann1.gbk','ann2.gb','ann3.gbf']:
        reffile = refd.joinpath(f)
        reffile.touch()
        if f.startswith('ref'):
            fastapaths.append(str(reffile.resolve()))
        elif f.startswith('ann'):
            annpaths.append(str(reffile.resolve()))
    seq_fofn.write_text('\n'.join(fastapaths))
    gbk_fofn.write_text('\n'.join(annpaths))
    return refd


@pytest.mark.parametrize("test_input,file_format", [
    (['.'],'fasta'),
    (['.'],'genbank'),
    (['ref1.fasta', 'ref2.fna', 'ref3.fa'],'fasta'),
    (['ann1.gbk', 'ann2.gb', 'ann3.gbf'], 'genbank'),
    (['refs.fofn'],'fasta'),
    (['anns.fofn'],'genbank'),
])
def test_file_list(test_input, file_format, refdir):
    expected = dict(
        fasta=['ref1.fasta','ref2.fna','ref3.fa'],
        genbank=['ann1.gbk','ann2.gb','ann3.gbf'],
    )
    os.chdir(refdir)
    assert sorted(fileManager.file_list(test_input, file_type=file_format)) == \
        [str(refdir.joinpath(_)) for _ in expected[file_format]]
