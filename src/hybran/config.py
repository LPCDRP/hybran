import tempfile
from types import SimpleNamespace

from Bio.Data import CodonTable


# TODO: update rest of code base to use this data
#       structure for tmpdir and other things,
#       like the generic ORF prefix.
cnf = SimpleNamespace(
    tmpdir=None,
    genetic_code=None,
    onegene=SimpleNamespace(
        min_identity=99,
        min_coverage=99,
    ),
    blast=SimpleNamespace(
        min_identity=95,
        min_coverage=95,
    ),
    mcl_inflation=1.5,
)

def init():
    # TODO: stop using this global variable and move towards using cnf.tmpdir
    global hybran_tmp_dir
    # Make the temporary files directory
    hybran_tmp_dir = tempfile.mkdtemp()
    cnf.tmpdir = hybran_tmp_dir

def set_genetic_code(code):
    cnf.genetic_code = code

def ratt(path, genetic_code, splice, correct_splice, correct_pseudo):
    """
    Write a RATT configuration file to `path` and using
    start and stop codons as in the genetic_code table

    :param path: path to filename in which to write config
    :param genetic_code: int genetic code table number
    :param splice: str splice donor/acceptor sequences. like XX..XX
    :param correct_splice: boolean
    :param correct_pseudo: boolean
    """

    switches = {
        True : '1',
        False : '0',
    }
    gcode = CodonTable.generic_by_id[genetic_code]
    with open(path, 'w') as ratt_config:
        ratt_config.write('#START\n')
        for codon in gcode.start_codons:
            ratt_config.write(codon + '\n')

        ratt_config.write('#STOP\n')
        for codon in gcode.stop_codons:
            ratt_config.write(codon + '\n')

        ratt_config.write('\n'.join(['#SPLICE',splice,'']))
        ratt_config.write('#CORRECTSPLICESITE\n')
        ratt_config.write(switches[correct_splice]+'\n')
        ratt_config.write('#CORRECTPSEUDOGENE\n')
        ratt_config.write(switches[correct_pseudo]+'\n')
