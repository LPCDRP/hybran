from copy import deepcopy
import functools
import logging
import os
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing

from . import (
    BLAST,
    config,
    designator,
)
from .annomerge import (
    fissionfuser,
    key_ref_gene,
    liftover_annotation,
)
from .bio import AutarkicSeqFeature, SeqIO
from .converter import convert_gbk_to_gff
from .demarcate import coord_check
from .lumberjack import (
    log_feature_fates,
    log_coord_corrections,
    log_pseudos,
)
from .pseudoscan import pseudoscan
from .util import mpbreakpoint


def postprocess(
        isolate_id,
        contigs,
        prokka_outdir,
        postprocess_outdir,
        ref_annotation,
        ref_proteome,
        seq_ident,
        seq_covg,
        nproc=1,
):

    logger = logging.getLogger('PostprocessAbInitio')

    invalid_features = []
    invalid_features_logfile = os.path.join(
        postprocess_outdir,
        'invalid_features.tsv',
    )
    corrected_orf_logfile = os.path.join(
        postprocess_outdir,
        'coord_corrections.tsv',
    )
    pseudoscan_logfile = os.path.join(
        postprocess_outdir,
        'pseudoscan_report.tsv',
    )
    os.makedirs(
        postprocess_outdir,
        exist_ok=True,
    )
    input_prokka_genbank = os.path.join(
        prokka_outdir,
        isolate_id + '.gbk',
    )

    prokka_features = {}
    prokka_records = []
    for i, record in enumerate(SeqIO.parse(input_prokka_genbank, 'genbank')):
        seqname  = '.'.join([isolate_id, contigs[i]])
        prokka_records.append(record)

        invalid_contig_features = postprocess_contig(
            seqname=seqname,
            record=record,
            ref_annotation=ref_annotation,
            ref_proteome=ref_proteome,
            seq_ident=seq_ident,
            seq_covg=seq_covg,
            nproc=nproc,
        )
        prokka_features[contigs[i]] = record.features
        invalid_features += invalid_contig_features
    out_gbk = os.path.join(
        postprocess_outdir,
        os.path.basename(input_prokka_genbank)
    )
    SeqIO.write(
        prokka_records,
        out_gbk,
        "genbank",
    )
    convert_gbk_to_gff(out_gbk)

    with open(invalid_features_logfile, 'w') as rejects_log:
        log_feature_fates(invalid_features, rejects_log)

    with open(corrected_orf_logfile, 'w') as corr_log:
        log_coord_corrections(prokka_features, corr_log)

    with open(pseudoscan_logfile, 'w') as p_log:
        log_pseudos(prokka_features, p_log)

    return prokka_features

def postprocess_contig(
        seqname,
        record,
        ref_annotation,
        ref_proteome,
        seq_ident,
        seq_covg,
        nproc=1,
):

    logger = logging.getLogger('PostprocessAbInitio')

    if not record.features:
        logger.warning(f"{seqname} HAS NO AB INITIO ANNOTATION")
        return []

    contig_features = []
    for f in record.features:
        f = AutarkicSeqFeature.fromSeqFeature(f)
        for part in f.location.parts:
            part.ref = seqname
        f.references = {seqname: record.seq}
        contig_features.append(f)

    logger.info(f'{seqname}: postprocessing ab initio CDS annotations using {nproc} process(es)')
    #
    # can contain results for hits to multiple reference genes
    abinit_blast_results_complete = {}
    # only contains results for the accepted reference gene hit
    abinit_blast_results = {}
    mp_postprocess_feature = functools.partial(
        postprocess_feature,
        ref_annotation=ref_annotation,
        ref_proteome=ref_proteome,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
    )
    prokka_contig_cdss = [f for f in contig_features if f.type == 'CDS']
    with multiprocessing.Pool(processes=nproc) as pool:
        contig_features, ref_matched, coords_corrected = zip(*pool.map(
            mp_postprocess_feature,
            contig_features,
        ))

    logger.info(f'{seqname}: {sum(ref_matched)} out of {len(prokka_contig_cdss)} ORFs matched to a reference gene')
    logger.info(f'{seqname}: Corrected coordinates for {sum(coords_corrected)} ab initio ORFs')

    logger.info(f"{seqname}: Checking for fragmented ab initio annotations")
    abinit_features_postprocessed_list, dropped_abinit_fragments = fissionfuser(
        contig_features,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
    )

    record.features = abinit_features_postprocessed_list

    logger.info(f"{seqname}: {len(dropped_abinit_fragments)} gene fragment pairs merged")
    logger.info(f"{seqname}: {len(abinit_features_postprocessed_list)} total remaining ab initio features")
    return dropped_abinit_fragments

def postprocess_feature(
        feature,
        ref_annotation,
        ref_proteome,
        seq_ident,
        seq_covg,
):
    """
    general postprocessing of individual ab initio features.
    Currently includes reference blastp matching and pseudoscan.

    :param feature: AutarkicSeqFeature ab initio feature
    :param ref_annotation: dict of curated reference annotations
    :param ref_proteome: str file name of multi fasta reference amino acid sequences
    :param seq_ident: int sequence identity percentage threshold for BLAST
    :param seq_covg: int sequence alignment coverage percent threshold for BLAST
    """
    ref_matched = False
    coords_corrected = False

    if 'gene' in feature.qualifiers:
        # When prokka assigns the same gene name to multiple orfs, it appends _1, _2, ... to make the names unique.
        # That causes issues for us because we expect all copies of a gene to have the same name.
        feature.qualifiers['gene'][0] = re.sub(r"_\d+$","",feature.qualifiers['gene'][0])

    if feature.type != 'CDS':
        return feature, ref_matched, coords_corrected

    top_hit, low_covg, blast_hits = BLAST.reference_match(
        SeqRecord(Seq(feature.qualifiers['translation'][0])),
        subject=ref_proteome,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
    )

    if top_hit:
        ref_matched = True
        ref_id, ref_ltag, ref_gene = top_hit.split('%%%')
        feature.source = ref_id
        feature.qualifiers['gene'] = [ref_gene]
        og_feature_location = deepcopy(feature.location)
        feature_is_pseudo = pseudoscan(
            feature,
            ref_annotation[key_ref_gene(ref_id, ref_gene)],
            seq_ident,
            seq_covg,
            attempt_rescue=True,
            blast_hit_dict=blast_hits[ref_gene]
        )

        coords_corrected = (og_feature_location != feature.location)

        liftover_annotation(
            feature,
            ref_annotation[key_ref_gene(ref_id, ref_gene)],
            inference=':'.join([
                f"similar to AA sequence",
                    ref_id,
                    ref_ltag,
                    ref_gene,
                    "blastp",
                ])
        )
        # TODO: if we still want to save blast results, we should make them attributes of the AutarkicSeqFeature object
        # abinit_blast_results[feature.qualifiers['locus_tag'][0]] = blast_hits[ref_gene]
    # Don't keep gene name assignments from Prokka. They can sometimes be based on
    # poor sequence similarity and partial matches (despite its --coverage option).
    # Keeping them is risky for propagation of the name during clustering.
    elif 'gene' in feature.qualifiers:
        designator.append_qualifier(
            feature.qualifiers, 'gene_synonym',
            feature.qualifiers['gene'][0],
        )
        feature.qualifiers.pop('gene', None)


    return feature, ref_matched, coords_corrected
