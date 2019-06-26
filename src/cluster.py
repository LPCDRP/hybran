#!/usr/bin/env python2.7
import argparse
import logging
import fastaFromGFF
import BLAST
import cdhit
import mcl


def arguments():
    parser = argparse.ArgumentParser(description='Clustering a set of sequences')
    parser.add_argument('-d', '--dir', help='Directory containing GFF annotations. Use independent of -f/--fasta')
    parser.add_argument('-f', '--fasta', help='CDS multiesequence FASTA to cluster. Use independent of -d/--dir')
    parser.add_argument('-n', '--nproc', help='Number of cores/processors to use. Default is 1',
                        default=1, type=int)
    parser.add_argument('-o', '--output', help='Output filename. Default is clustered_proteins',
                        default='clustered_proteins')
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s')
    logger = logging.info('__main__')
    args = arguments()
    if args.dir:
        gff_gene_dict = fastaFromGFF.create_fasta(directory=args.dir)
        fasta = 'cds_seqs.fasta'
    elif args.fasta:
        fasta = args.fasta
        gff_gene_dict = ''
    else:
        print 'Did not provide -d/--dir or -f/--fasta. A choice is required. Exiting...'
        exit()
    clusters = cdhit.cd_hit(nproc=args.nproc,
                            fasta=fasta,
                            out='cdhit_clusters.fasta')
    BLAST.run_blast(fastafile='cdhit_clusters.fasta')
    mcl.run_mcl(in_blast='blast_results',
                cdhit_clusters=clusters,
                out_name=args.output,
                gene_names=gff_gene_dict)


if __name__ == '__main__':
    main()
