import tempfile

from Bio.Data import CodonTable

def init():
    global hybran_tmp_dir
    # Make the temporary files directory
    hybran_tmp_dir = tempfile.mkdtemp()

def set_genetic_code(code):
    global genetic_code
    genetic_code = code

def ratt(path, genetic_code):
    """
    Write a RATT configuration file to `path` and using
    start and stop codons as in the genetic_code table

    :param path: path to filename in which to write config
    :param genetic_code: int genetic code table number
    """

    gcode = CodonTable.generic_by_id[genetic_code]
    with open(path, 'w') as ratt_config:
        ratt_config.write('#START\n')
        for codon in gcode.start_codons:
            ratt_config.write(codon + '\n')

        ratt_config.write('#STOP\n')
        for codon in gcode.stop_codons:
            ratt_config.write(codon + '\n')

        ratt_config.write('#SPLICE\nXX..XX\n')
        ratt_config.write('#CORRECTSPLICESITE\n')
        ratt_config.write('1\n')
        ratt_config.write('#CORRECTPSEUDOGENE\n')
        ratt_config.write('0\n')
