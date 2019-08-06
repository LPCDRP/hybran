#!/usr/bin/python2.7
import os
import argparse
import subprocess
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def run_ratt_prokka_annomerge(ref_dir, fasta, ref_cds, script_dir, cpus):
    c = os.getcwd()
    isolate = fasta.split('/')[-1].split('.')[0]
    try:
        os.mkdir(isolate)
    except OSError:
        pass
    os.chdir(isolate)
    cmd = [script_dir + 'ratt_prokka.sh',
           ref_dir,
           fasta,
           isolate,
           ref_cds,
           cpus]
    print ' '.join(cmd)
    subprocess.Popen(cmd)
    os.chdir(c)


def get_first_reference_proteome(genbank):
    logger = logging.getLogger('ReferenceDatabase')
    logger.info('Creating a reference proteome FASTA for Prokka from ' + genbank)
    seqs = []
    for record in SeqIO.parse(genbank, 'genbank'):
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS':
                    seq_record = SeqRecord(feature.qualifiers['translation'],
                                           id='1',
                                           description='~~~' + feature.qualifiers['gene'][0] + '~~~')
                    seqs.append(seq_record)
    with open('ref_proteome.fasta', 'w') as output:
        for s in seqs:
            SeqIO.write(s, output, 'fasta')
    return 'ref_proteome.fasta'


def cmds():
    parser = argparse.ArgumentParser(description='Annotate TUBerculosis: a pipeline to annotate Mycobacterium '
                                                 'tuberculosis de novo assembled genomes. Annotation of other species '
                                                 'or mixing species within an annotation run is NOT recommended.')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    required.add_argument('-g', '--genomes', help='Directory containing all genomes desired to be annotationed. '
                                                  'FASTA format required')
    required.add_argument('-r', '--references', help='Directory containing EMBL and Genbank files of reference '
                                                     'annotations to transfer. Only the first 30 reference annotations '
                                                     'will be transferred with RATT and the first annotation will '
                                                     'be used as the reference database in the Prokka reference '
                                                     'step.')
    optional.add_argument('-e', '--emapper', help='Directory of the eggnog-mapper repository. Default is the current '
                                                  'working directory. Full path only',
                          default=os.getcwd())
    optional.add_argument('-o', '--output', help='Directory to output all new annotation files. Default is the current '
                                                 'working directory. Full path only',
                          default=os.getcwd())
    optional.add_argument('-n', '--nproc', help='Number of processors/CPUs to use. Default is 1',
                          default='1')
    return parser.parse_args()


def main():
    args = cmds()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s')
    logger = logging.getLogger('AnnoTUB')
    # Verify the RATT_HOME variable is set
    if not os.environ['RATT_HOME']:
        logger.error('RATT_HOME variable is not defined. '
                     'Please export RATT_HOME=[location of downloaded ratt executable/config')
        exit(EnvironmentError)
    # Verify eggnog-mapper/ is present in given directory
    if not os.path.isdir(args.emapper):
        logger.error('eggnog-mapper is not in ' + args.emapper + '. Location of eggnog-mapper/ repo is required')
        exit(OSError)
    cwd = os.getcwd() + '/'
    if not args.genomes.startswith('/'):
        args.genomes = cwd + args.genomes
    if not args.references.startswith('/'):
        args.references = cwd + args.references
    if not args.output.startswith('/'):
        args.output = cwd + args.output
    os.chdir(args.output)
    refdir = args.output + 'temp_references/'
    try:
        os.mkdir(refdir)
    except OSError:
        pass
    logger.info('Getting first 30 reference annotations')
    embls = [i for i in os.listdir(args.references) if i.endswith('.embl')][0:29]
    for e in embls:
        try:
            os.symlink(args.references + e, refdir + e)
        except OSError:
            continue
    ref_cds = get_first_reference_proteome(args.references + embls[0].split('.')[0] + '.gbk')
    for f in os.listdir(args.genomes):
        if f.endswith('.fasta'):
            run_ratt_prokka_annomerge(ref_dir=refdir,
                                      fasta=args.genomes + f,
                                      ref_cds=ref_cds,
                                      script_dir=cwd,
                                      cpus=args.nproc)


if __name__ == '__main__':
    main()
