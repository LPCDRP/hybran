import functools
import os
import logging

import Bio
from Bio.Seq import Seq, translate
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord
# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing

from . import BLAST
from . import config
from . import converter
from . import designator
from . import extractor
from .annomerge import coord_check
from .annomerge import fusionfisher
from .annomerge import get_and_remove_ref_tracer
from .annomerge import get_ordered_features
from .annomerge import has_broken_stop
from .annomerge import key_ref_gene
from .annomerge import pseudoscan
from .bio import AutarkicSeqFeature, SeqIO
from .lumberjack import log_feature_fate
from .lumberjack import log_coord_corrections
from .lumberjack import log_pseudos
from .util import mpbreakpoint

def postprocess(
        isolate_id,
        contigs,
        ratt_outdir,
        postprocess_outdir,
        ref_annotation,
        seq_ident,
        seq_covg,
        nproc=1,
        enforce_thresholds=False,
):

    logger = logging.getLogger('PostprocessRATT')

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

    global genetic_code
    genetic_code = config.genetic_code

    ratt_correction_files = []
    ratt_features = {}
    try:
        for contig in contigs:
            seqname = '.'.join([isolate_id, contig])
            embl_file = f"{isolate_id.replace('|','_')}.{contig.replace('|','_')}.final.embl"
            gbk = converter.convert_embl_to_gbk(
                os.path.join(ratt_outdir, embl_file),
            )
            ratt_contig_record, invalid_contig_features = postprocess_contig(
                seqname=seqname,
                ratt_gbk=gbk,
                ref_annotation=ref_annotation,
                seq_ident=seq_ident,
                seq_covg=seq_covg,
                nproc=nproc,
                enforce_thresholds=enforce_thresholds,
            )
            ratt_features[contig] = ratt_contig_record.features
            SeqIO.write(
                ratt_contig_record,
                os.path.join(postprocess_outdir, os.path.basename(gbk)),
                "genbank",
            )
            invalid_features += invalid_contig_features

    except OSError:
        logger.error('Expecting RATT annotation files but found none')
    if not ratt_features:
        logger.error('RATT did not complete running. Please see the log for more details.')


    with open(invalid_features_logfile, 'w') as rejects_log:
        [log_feature_fate(_[0], rejects_log, _[1]) for _ in invalid_features]

    with open(corrected_orf_logfile, 'w') as corr_log:
        log_coord_corrections(ratt_features, corr_log)

    with open(pseudoscan_logfile, 'w') as p_log:
        log_pseudos(ratt_features, p_log)

    return ratt_features

def postprocess_contig(
        seqname,
        ratt_gbk,
        ref_annotation,
        seq_ident,
        seq_covg,
        nproc=1,
        enforce_thresholds=False,
):

    logger = logging.getLogger('PostprocessRATT')

    ratt_contig_record = SeqIO.read(ratt_gbk, 'genbank')

    ratt_contig_features = []
    for feature in ratt_contig_record.features:
        feature = AutarkicSeqFeature.fromSeqFeature(feature)
        ref_contig_id = get_and_remove_ref_tracer(feature)
        feature.source = ref_contig_id
        if feature.location is not None:
            for part in feature.location.parts:
                part.ref = seqname
        feature.references = {seqname: ratt_contig_record.seq}
        # maybe RATT should be adding this inference tag itself
        if 'locus_tag' in feature.qualifiers:
            infer_string = ':'.join([
                f"similar to nucleotide sequence",
                ref_contig_id,
                extractor.get_ltag(feature),
                extractor.get_gene(feature),
                "RATT",
            ])
        else:
            infer_string = ':'.join([
                f"similar to nucleotide sequence",
                ref_contig_id,
                "RATT",
            ])
        designator.append_qualifier(feature.qualifiers, 'inference', infer_string)
        if feature.type not in [
                'source',
                'gene',
        ]:
            ratt_contig_features.append(feature)


    logger.info(f"Validating CDSs using {nproc} process(es)")
    valid_features = []
    invalid_ratt_features = []
    mp_validate = functools.partial(
        validate,
        ref_annotation=ref_annotation,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
        enforce_thresholds=enforce_thresholds,
    )
    with multiprocessing.Pool(processes=nproc) as pool:
        results = pool.map(
            mp_validate,
            ratt_contig_features,
        )
    for valid, feature, remark in results:
        valid_features.append(feature) if valid else invalid_ratt_features.append((feature, remark))

    logger.info(f"{seqname}: {len(invalid_ratt_features)} RATT features failed validation.")
    ratt_contig_features = get_ordered_features(valid_features)

    logger.info(f"{seqname}: Checking for gene fusion signatures in RATT annotations...")
    ratt_contig_features, merged_features, inconsistent_ratt_features = fusionfisher(ratt_contig_features)
    logger.info(f"{seqname}: {len(merged_features)} RATT fusion genes detected.")
    logger.info(f"{seqname}: {len(inconsistent_ratt_features)} RATT annotations found to be inconsistent.")

    ratt_rejects = invalid_ratt_features + inconsistent_ratt_features
    ratt_contig_record.features = ratt_contig_features

    logger.info(f"{seqname}: {len(ratt_contig_features)} total remaining RATT features")

    # if merged_features:
    #     merged_features_record = ratt_contig_record[:]
    #     merged_features_record.features = merged_features
    #     SeqIO.write(merged_features_record, output_merged_genes, 'genbank')

    return ratt_contig_record, ratt_rejects

def validate(
        feature,
        ref_annotation,
        seq_ident,
        seq_covg,
        enforce_thresholds,
):
    """
    This function takes as input a list of features and checks if the length of the CDSs are divisible by
    3 and if the CDS is split across multiple locations. If so, it outputs the features to stdout and removes them
    from the valid_ratt_annotations. The function BLASTs the sequence to corresponding amino acid sequence in
    reference as well.
    :param feature: AutarkicSeqFeature feature annotation from RATT
    :return: List of valid RATT features (list of SeqFeature objects)
    """
    logger = logging.getLogger('ValidateRATTCDSs')
    blast_stats = {}
    valid = False
    unbroken = False
    remark = ''

    if enforce_thresholds:
        ratt_seq_ident = seq_ident
        ratt_seq_covg = seq_covg
    else:
        ratt_seq_ident = ratt_seq_covg = 0

    if feature.location is None:
        valid = False
        remark = "Empty feature location"
        return valid, feature, remark

    if feature.type != 'CDS':
        if feature.type in [
                'rRNA', # these aren't reliably transferred by RATT
                'tRNA',
        ]:
            valid = False
            remark = f"{feature.type}s categorically rejected in favor of ab initio"
        else:
            valid = True
        return valid, feature, remark

    compound_interval = isinstance(feature.location,Bio.SeqFeature.CompoundLocation)
    # Identify features with 'joins'
    if compound_interval and 'ribosomal_slippage' not in feature.qualifiers:
        #Need to initialize the feature without the compound location attribute.
        #The earliest start and the latest end of the joined feature will be bridged together
        feature_start = feature.location.start
        feature_end = feature.location.end
        feature_strand = feature.strand
        feature.location = SimpleLocation(
            feature_start,
            feature_end,
            strand=feature_strand,
            ref=feature.location.parts[0].ref,
        )
        #Check if feature has an internal stop codon.
        #
        # If it doesn't, we will assign pseudo and accept it.

        # The gff conversion of a gbk entry with joins is not meaningful,
        # and causes some problems, as the entire sequence gets labeled
        # "biological region" and two basically empty CDS records are created.

        broken_stop, stop_note = has_broken_stop(feature)
        if broken_stop:
            good_start, good_stop = coord_check(
                feature,
                ref_annotation[key_ref_gene(feature.source, feature.qualifiers['gene'][0])],
                fix_start=True,
                fix_stop=True,
            )
            if not good_stop:
                remark = "RATT-introduced compound interval did not include reference stop position."
                valid = False
                return valid, feature, remark

    feature_is_pseudo = pseudoscan(
        feature,
        ref_annotation[key_ref_gene(feature.source, feature.qualifiers['gene'][0])],
        seq_ident,
        seq_covg,
        attempt_rescue=True,
    )

    if feature_is_pseudo:
        valid = True
    else:
        unbroken = True

    if unbroken:
        feature_sequence = translate(
            feature.extract(),
            table=genetic_code,
            to_stop=True,
        )
        # TODO: if we want to keep track of the blast stats, we should add it as an attribute to the AutarkicSeqFeature object
        #ratt_blast_results.update(blast_stats)
        feature.qualifiers['translation'] = [str(feature_sequence)]

        if not enforce_thresholds:
            valid = True
        else:
            ref_feature = ref_annotation[
                key_ref_gene(feature.source, feature.qualifiers['gene'][0])
            ]
            try:
                ref_seq = ref_feature.qualifiers['translation'][0]
            except KeyError:
                ref_seq = translate(
                    extractor.get_seq(ref_feature),
                    table=genetic_code, to_stop=True
                )

            if len(feature_sequence) == 0:
                remark = 'length of AA sequence is 0'
            else:
                top_hit, low_covg, blast_stats = BLAST.reference_match(
                    query=SeqRecord(feature_sequence),
                    subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                    seq_ident=seq_ident,
                    seq_covg=seq_covg,
                )

                if top_hit:
                    valid = True
                else:
                    remark = 'No blastp hit to corresponding reference CDS at specified thresholds.'


    return valid, feature, remark
