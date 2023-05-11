import logging
import os
import subprocess
from urllib.error import HTTPError

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

from . import designator

def get_ltag(feature):
    return feature.qualifiers['locus_tag'][0]

def get_gene(feature):
    try:
        gene = feature.qualifiers['gene'][0]
    except KeyError:
        gene = get_ltag(feature)
    return gene

def get_seq(feature):
    """
    Adapted from Bio.SeqFeature CompoundLocation.extract()
    """
    parts = [
        loc.extract(
            feature.references[loc.ref],
            references=feature.references
        ) for loc in feature.location.parts
    ]
    return functools.reduce(lambda x, y: x + y, parts)

def get_contig_names(fasta_file):
    return [record.id for record in SeqIO.parse(fasta_file,"fasta")]

def get_genetic_code(genbank):
    """
    Find genetic code from annotation file
    """

    # default to standard bacterial genetic code
    # if we don't find anything different.
    gcode = 11

    for record in SeqIO.parse(genbank, 'genbank'):
        if record.features:
            for feature in record.features:
                if 'transl_table' in feature.qualifiers:
                    gcode = int(feature.qualifiers['transl_table'][0])
                    break

    return gcode

def gene_dict(genbank, cds_only=True):
    """
    Create a dictionary mapping locus tags to gene names

    :param genbank: genbank annotation file name
    :param cds_only: bool whether to only consider coding sequences
    :return: dict
    """

    genes = dict()

    for record in SeqIO.parse(genbank, 'genbank'):
        if record.features:
            for feature in record.features:
                if((not cds_only and 'locus_tag' in feature.qualifiers.keys())
                   or (cds_only and feature.type == 'CDS')
                ):
                    genes[get_ltag(feature)] = get_gene(feature)
    return genes

def get_taxonomy_id(genbank):
    """
    Attempt to find the NCBI Taxonomy ID from the genbank accession number

    :param genbank: reference genbank annotation file
    :return: str tax_id corresponding to the annotation file, or None
             if the query failed.
    """
    logger = logging.getLogger('GetReferenceTaxonomyID')
    gb = list(SeqIO.parse(genbank, format="genbank"))[0]
    Entrez.email = 'A.N.Other@example.com'
    tax_id = None

    # if organism isn't specific (strain level), this won't be what we want
    #organism = gb.description
    #tax_id = Entrez.read(
    #    Entrez.esearch(db = "Taxonomy", term = organism)
    #)

    accession = gb.id
    try:
        nuc_record = Entrez.read(
            Entrez.efetch(
                db="nucleotide",
                id=accession,
                retmode="xml",
            )
        )[0]
        tax = [_['GBQualifier_value'] for _ in nuc_record['GBSeq_feature-table'][0]['GBFeature_quals']
               if _['GBQualifier_value'].startswith('taxon:')][0]
        if tax:
            tax_id = tax.split(':')[1]
    except HTTPError:
        pass

    if tax_id:
        logger.info("Using taxonomy ID " + str(tax_id) + " for " + accession)
    else:
        logger.warning("Could not find a taxonomy ID for '" + accession + "'.")

    return tax_id

def prokka_faa(feature):
    """
    Create fasta descriptions in Prokka format.
    See https://github.com/tseemann/prokka#fasta-database-format
    :param feature: sequence feature
    :return: string representing fasta record description in prokka format
    """
    # Don't include EC for now.
    # Not sure how to handle multiple EC numbers in Prokka faa format
    # ec = feature.qualifiers['EC_number'][0]
    ec = ''
    try:
        product = feature.qualifiers['product'][0]
    except KeyError:
        product = ''
    return '~~~'.join([ec,get_gene(feature),product])

def fastaFromGbk(genbank, out_cds, out_genome,
                 identify = lambda f: ':'.join([get_ltag(f),
                                                get_gene(f)]),
                 describe = lambda f: '',
                 ):
    """
    Extracts amino acid CDS sequences from a Genbank annotation file

    :param genbank:
    :param out_cds: file name or handle in which to write CDS sequences
    :param out_genome: file name or handle in which to write genome sequence
    :param identify: function to apply to feature record to get fasta record ID
    :param describe: function to apply to feature record to get fasta record description
    :return:
    """
    logger = logging.getLogger('FastaFromGbk')
    contigs = []
    seqs = []
    n_named_cds = 0
    ref_id = os.path.basename(genbank)
    # this is a bit circular, but I don't want to deal with the first CDS
    # possibly being a pseudogene.
    genetic_code = get_genetic_code(genbank)
    for record in SeqIO.parse(genbank, 'genbank'):
        ref_contig_id = '.'.join([ref_id, record.name])
        contigs.append(SeqRecord(
            record.seq,
            id=record.id, # os.path.splitext(os.path.basename(genbank))[0],
            description=''
        ))
        if record.features:
            for feature in record.features:
                feature.ref = ref_contig_id
                if feature.type == 'CDS':
                    if (
                            'gene' in feature.qualifiers
                            and not designator.is_unannotated(feature.qualifiers['gene'][0])
                            and feature.qualifiers['gene'][0] != feature.qualifiers['locus_tag'][0]
                    ):
                        n_named_cds += 1
                    if designator.is_pseudo(feature.qualifiers):
                        seq_record = SeqRecord(
                            translate(feature.extract(record.seq), table=genetic_code, to_stop=True),
                            id=identify(feature),
                            description=describe(feature))
                    else:
                        seq_record = SeqRecord(
                            Seq(feature.qualifiers['translation'][0]),
                            id=identify(feature),
                            description=describe(feature))
                    seqs.append(seq_record)
    SeqIO.write(seqs, out_cds, 'fasta')
    SeqIO.write(contigs, out_genome, 'fasta')
    return n_named_cds

def subset_fasta(inseq, outseq, match, identify = lambda _:_):
    """
    write a new fasta file containing only sequences with
    matching record IDs.
    :param infile: str input fasta file name
    :param outseq: str output fasta file name
    :param match: function to apply to record.id for a Boolean result
    :param identify: function to transform record.id prior to matching
    """
    seqs = []
    for record in SeqIO.parse(inseq, 'fasta'):
        record.id = identify(record.id)
        if match(record.id):
            seqs.append(record)
    SeqIO.write(seqs, outseq, 'fasta')

def grep_seqs(gff):
    gff_lines = []
    cmd = ['grep', 'translation=', gff]
    translations = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in translations.stdout:
        try:
            gff_lines.append(line)
        except TypeError:
            continue
    return gff_lines


def fastaFromGffList(gffs, out_cds):
    """
    Extracts amino acid CDS sequences from GFF annotation files

    :param gffs: list of GFF annotation file names
    :param out_cds: file name or handle in which to write CDS sequences
    :return gff_gene_dict: dictionary of CDS sequences by sample-recordID
    """
    logger = logging.getLogger('FastaFromGff')
    seqs = []
    gff_gene_dict = {}
    for gff in gffs:
        raw_out = grep_seqs(gff)
        gff_name = os.path.splitext(os.path.basename(gff))[0]
        for line in raw_out:
            gff_id = gff_name + '|' + [j.split('=')[1] for i in line.split('\t')
                                       if i.startswith('ID=') for j in i.split(';')][0]
            gene = None
            for i in line.split(';'):
                if i.startswith('gene='):
                    gene = [i.split('=')[1].rstrip('\n')][0]
                if i.startswith('locus_tag='):
                    locus_tag = [i.split('=')[1].rstrip('\n')][0]
            if not gene:
                gene = locus_tag
            gff_gene_dict[gff_id] = gene
            translation = [i.split('=')[1] for i in line.split(';') if i.startswith('translation=')][0]
            record = SeqRecord(Seq(translation.rstrip('\n')),
                               id=gff_id,
                               description='')
            seqs.append(record)

    SeqIO.write(seqs, out_cds, 'fasta')
    return gff_gene_dict
