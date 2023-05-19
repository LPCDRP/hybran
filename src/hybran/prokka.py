from copy import deepcopy
import functools
import logging
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing

from . import BLAST
from . import config
from . import designator
from .annomerge import coord_check
from .annomerge import fissionfuser
from .annomerge import key_ref_gene
from .annomerge import liftover_annotation
from .annomerge import log_feature_fate
from .annomerge import pseudoscan
from .bio import AutarkicSeqFeature


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
    os.makedirs(
        postprocess_outdir,
        exist_ok=True,
    )
    input_prokka_genbank = os.path.join(
        prokka_outdir,
        isolate_id + '.gbk',
    )

    global genetic_code
    genetic_code = config.genetic_code

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
    SeqIO.write(
        record,
        os.path.join(
            postprocess_outdir,
            os.path.basename(input_prokka_genbank)
        ),
        "genbank",
    )

    with open(invalid_features_logfile, 'w') as rejects_log:
        [log_feature_fate(_[0], rejects_log, _[1]) for _ in invalid_features]

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

    contig_features = []
    for f in record.features:
        f = AutarkicSeqFeature.fromSeqFeature(f)
        for part in f.location.parts:
            part.ref = seqname
        f.references = {seqname: record.seq}
        if 'gene' in f.qualifiers.keys():
            # When prokka assigns the same gene name to multiple orfs, it appends _1, _2, ... to make the names unique.
            # That causes issues for us because we expect all copies of a gene to have the same name.
            f.qualifiers['gene'][0] = re.sub(r"_\d+$","",f.qualifiers['gene'][0])
        contig_features.append(f)


    record.features = contig_features


    logger.info(f'{seqname}: Checking ab initio CDS annotations for matches to reference using {nproc} process(es)')
    #
    # can contain results for hits to multiple reference genes
    abinit_blast_results_complete = {}
    # only contains results for the accepted reference gene hit
    abinit_blast_results = {}
    refmatch = functools.partial(
        BLAST.reference_match,
        subject=ref_proteome,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
    )
    prokka_contig_cdss = [f for f in contig_features if f.type == 'CDS']
    with multiprocessing.Pool(processes=nproc) as pool:
        blast_package = pool.map(
            refmatch,
            [SeqRecord(Seq(f.qualifiers['translation'][0])) for f in prokka_contig_cdss],
        )
    n_coords_corrected = 0
    n_bad_starts = 0
    for j in range(len(prokka_contig_cdss)):
        top_hit, low_covg, blast_hits = blast_package[j]
        feature = prokka_contig_cdss[j]
        if top_hit:
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

            if (og_feature_location != feature.location):
                n_coords_corrected += 1

                # # for logging purposes
                # if feature_is_pseudo:
                #     corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']
                #     corrected_orf_report[-1][1].qualifiers['pseudo'] = ['']
                # else:
                #     corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']

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
            abinit_blast_results[feature.qualifiers['locus_tag'][0]] = blast_hits[ref_gene]
        # Don't keep gene name assignments from Prokka. They can sometimes be based on
        # poor sequence similarity and partial matches (despite its --coverage option).
        # Keeping them is risky for propagation of the name during clustering.
        elif 'gene' in feature.qualifiers:
            designator.append_qualifier(
                feature.qualifiers, 'gene_synonym',
                feature.qualifiers['gene'][0],
            )
            feature.qualifiers.pop('gene', None)
        # We save all the blast hits at this point in case coordinates were corrected and the hits changed
        abinit_blast_results_complete[feature.qualifiers['locus_tag'][0]] = blast_hits

    logger.debug(f'{seqname}: {len(abinit_blast_results.keys())} out of {len(prokka_contig_cdss)} ORFs matched to a reference gene')
    logger.debug(f'{seqname}: Corrected coordinates for {n_coords_corrected} ab initio ORFs')

    logger.info(f"{seqname}: Checking for fragmented ab initio annotations")
    abinit_features_postprocessed_list, dropped_abinit_fragments = fissionfuser(
        contig_features,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
        abinit_blast_results=abinit_blast_results,
    )
    logger.debug(f"{seqname}: {len(dropped_abinit_fragments)} gene fragment pairs merged")
    return dropped_abinit_fragments
