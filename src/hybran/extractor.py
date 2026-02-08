import functools
import logging
import os
import re
import subprocess
from urllib.error import HTTPError

import Bio
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from . import designator
from .bio import (
    SeqIO,
    translate,
)
from .util import keydefaultdict


def get_ltag(feature):
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    else:
        return None

def get_gene(feature, tryhard=True):
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif tryhard:
        return get_ltag(feature)
    else:
        return None

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

def get_source(feature):
    """
    Get the reference genome identifier for a (Autarkic)SeqFeature.
    If the `source` attribute is not set, assume the feature itself is the reference and return its own record ID.
    """
    if hasattr(feature, 'source') and feature.source:
        return feature.source
    else:
        return feature.location.parts[0].ref

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

def gene_dict(features_by_contig, cds_only=True):
    """
    Create a dictionary mapping locus tags to gene names

    :param genbank: genbank annotation file name
    :param cds_only: bool whether to only consider coding sequences
    :return: dict
    """

    genes = dict()

    for record in features_by_contig:
        for feature in features_by_contig[record]:
            if (
                    (not cds_only and 'locus_tag' in feature.qualifiers)
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
                 identify = lambda f: '%%%'.join([get_ltag(f),
                                                get_gene(f)]),
                 describe = lambda f: '',
                 ):
    """
    Extracts amino acid CDS sequences from a Genbank annotation file

    :param genbank:
    :param out_cds:
        file name or handle in which to write CDS sequences
        If None, CDS file will not be written.
    :param out_genome:
        file name or handle in which to write genome sequence
        If None, genome sequence file will not be written.
    :param identify: function to apply to feature record to get fasta record ID
    :param describe: function to apply to feature record to get fasta record description
    :return:
        contigs     - list of SeqRecords for contig sequences
        cds_seqs    - list of SeqRecords for CDS amino acid sequences
        n_named_cds - int number of named CDSs
    """
    logger = logging.getLogger('FastaFromGbk')
    contigs = []
    seqs = []
    n_named_cds = 0
    ref_id = os.path.basename(os.path.splitext(genbank)[0])
    # this is a bit circular, but I don't want to deal with the first CDS
    # possibly being a pseudogene.
    genetic_code = get_genetic_code(genbank)
    for record in SeqIO.parse(genbank, 'genbank'):
        ref_contig_id = '.'.join([ref_id, record.id])
        contigs.append(SeqRecord(
            record.seq,
            id=record.id, # os.path.splitext(os.path.basename(genbank))[0],
            description=''
        ))
        if record.features:
            for feature in record.features:
                for part in feature.location.parts:
                    part.ref = ref_contig_id
                feature.references = {ref_contig_id: record.seq}
                if feature.type == 'CDS':
                    if (
                            'gene' in feature.qualifiers
                            and not designator.is_unannotated(feature.qualifiers['gene'][0])
                            and feature.qualifiers['gene'][0] != feature.qualifiers['locus_tag'][0]
                    ):
                        n_named_cds += 1
                    if 'translation' not in feature.qualifiers:
                        seq_record = SeqRecord(
                            translate(get_seq(feature), table=genetic_code, to_stop=True),
                            id=identify(feature),
                            description=describe(feature))
                    else:
                        seq_record = SeqRecord(
                            Seq(feature.qualifiers['translation'][0]),
                            id=identify(feature),
                            description=describe(feature))
                    seqs.append(seq_record)
    if out_cds:
        SeqIO.write(seqs, out_cds, 'fasta')
    if out_genome:
        SeqIO.write(contigs, out_genome, 'fasta')
    return contigs, seqs, n_named_cds

def subset_fasta(inseq, outseq, match, identify = lambda _:_):
    """
    write a new fasta file containing only sequences with
    matching record IDs.
    :param infile: str input fasta file name
    :param outseq: str output fasta file name
    :param match: function to apply to record.id for a Boolean result
    :param identify: function to transform record.id prior to matching
    :return: set unique_ids (matched and post-transformation with identify())
    """
    seqs = []
    unique_ids = set()
    for record in SeqIO.parse(inseq, 'fasta'):
        record.id = identify(record.id)
        if match(record.id):
            seqs.append(record)
            unique_ids.add(record.id)
    if outseq:
        SeqIO.write(seqs, outseq, 'fasta')
    return unique_ids

def grep_seqs(gff):
    """
    Get all GFF annotations that contain a translation

    :param gff: str GFF file name
    :return: list of str lines from the gff file that contain translations.
    """
    gff_lines = []
    cmd = ['grep', 'translation=', gff]
    translations = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in translations.stdout:
        try:
            gff_lines.append(line)
        except TypeError:
            continue
    return gff_lines

def get_and_remove_ref_tracer(feature):
    """
    Remove the temporary tracer note we added to keep track of where an annotation that RATT placed originated from
    """
    ref_contig_id = ""
    if 'note' in feature.qualifiers:
        for note in feature.qualifiers['note']:
            if note.startswith("HYBRANSOURCE"):
                ref_contig_id = re.sub(r'\s+', '', note.split(':')[1])
                feature.qualifiers['note'].remove(note)
                break

    return ref_contig_id

def ref_fuse(fusion_gene_name, annotations):
    """
    Create a dummy SeqFeature to use as a reference for fusion gene coordinate checking/correction.
    Throws a KeyError if one of the constituent names is not in the annotations dictionary or if the fusion gene name is not valid.
    This is intended to be used as a callable for defaultdict for reference annotations.

    :param fusion_gene_name: A string, expected to be @@@-delimited between ref and gene, and each of those to be ::-delimited, corresponding to the ref and gene name of each component.
    :return: SeqFeature with concatenated coordinates
    """

    (refs, fusion_name) = fusion_gene_name.split('@@@')
    const_genes = fusion_name.split('::')
    const_refs = refs.split('::')
    # Not a fusion gene; defaultdict lookup should fail
    if len(const_genes) <= 1:
        raise KeyError(f'no gene "{const_genes[0]}" found for reference "{const_refs[0]}"') from None
    location_parts = []
    location_sequences = {}
    for i in range(len(const_genes)):
        ref_feature_i = annotations[designator.key_ref_gene(const_refs[i], const_genes[i])]
        location_sequences.update(ref_feature_i.references)
        location_parts += list(ref_feature_i.location.parts)
    ref_fusion = SeqFeature(
        # Biopython currently doesn't support the 'order' operator for feature.extract()
        Bio.SeqFeature.CompoundLocation(location_parts, operator='join'),
        qualifiers={'locus_tag':fusion_gene_name, 'gene':fusion_gene_name},
    )

    ref_fusion.references = location_sequences
    return ref_fusion

def load_gbks(gbk_list, feature_types=['CDS']):
    feature_dict = keydefaultdict(ref_fuse)
    for gbk in gbk_list:
        feature_dict.update(load_gbk(gbk, feature_types=feature_types))
    return feature_dict

def load_gbk(gbk_file, feature_types=['CDS']):
    """
    Load a reference annotation (after some light postprocessing) into a flat dictionary indexed by the reference contig ID joined with the gene name.
    Note that a gene name is assumed to be present for all features, an assumption which is fulfilled by our reference annotation preprocessing.

    :param gbk_file: str path to genbank annotation file to load
    :param feature_types: list of annotation types to include.
    """
    feature_dict = keydefaultdict(ref_fuse)
    strain_id = os.path.splitext(os.path.basename(gbk_file))[0]
    for record in SeqIO.parse(gbk_file, "genbank"):
        contig_id = '.'.join([strain_id, record.id])
        for feature in record.features:
            get_and_remove_ref_tracer(feature) # prevent our tracer note from propagating to future liftovers
            if feature.type not in feature_types:
                continue
            # setting feature.ref doesn't work for CompoundLocations
            for part in feature.location.parts:
                part.ref = contig_id
            feature.references = {contig_id: record.seq}
            # if reference paralogs have been collapsed, the last occurrence in the genome
            # will prevail.
            feature_dict[designator.key_ref_gene(contig_id, feature.qualifiers['gene'][0])] = feature
    return feature_dict

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
            line = line.rstrip()
            gff_id = gff_name + '@@@' + [j.split('=')[1] for i in line.split('\t')
                                       if i.startswith('ID=') for j in i.split(';')][0]
            gene = None
            for i in line.split(';'):
                if i.startswith('gene='):
                    gene = i.split('=')[1]
                elif i.startswith('locus_tag='):
                    locus_tag = i.split('=')[1]
                elif i.startswith('translation='):
                    translation = i.split('=')[1]
            # intermediate locus tags will have been copied
            # to the gene name field by this point and final
            # locus tags (not useful for determining the annotation
            # status) will have already been set.
            if not gene:
                continue
            gff_gene_dict[gff_id] = gene
            record = SeqRecord(Seq(translation),
                               id=gff_id,
                               description='')
            seqs.append(record)

    SeqIO.write(seqs, out_cds, 'fasta')
    return gff_gene_dict
