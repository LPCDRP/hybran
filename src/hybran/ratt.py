import functools
import os
import logging

import Bio
from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord
# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing

from . import (
    BLAST,
    config,
    converter,
    designator,
    extractor,
)
from .annomerge import (
    fusionfisher,
    get_and_remove_ref_tracer,
    get_ordered_features,
    key_ref_gene,
)
from .bio import (
    AutarkicSeqFeature,
    FeatureProperties,
    SeqIO,
    translate,
)
from .config import cnf
from .demarcate import coord_check, has_broken_stop
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
            out_gbk = os.path.join(postprocess_outdir, os.path.basename(gbk))
            SeqIO.write(
                ratt_contig_record,
                out_gbk,
                "genbank",
            )
            converter.convert_gbk_to_gff(out_gbk)
            invalid_features += invalid_contig_features

    except OSError:
        logger.error('Expecting RATT annotation files but found none')
    if not ratt_features:
        logger.error('RATT did not complete running. Please see the log for more details.')


    with open(invalid_features_logfile, 'w') as rejects_log:
        log_feature_fates(invalid_features, rejects_log)

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

    if not ratt_contig_record.features:
        logger.warning(f"{seqname} HAS NO RATT ANNOTATION")
        return ratt_contig_record, []

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
    for valid, feature, evid, remark in results:
        valid_features.append(feature) if valid else invalid_ratt_features.append({
            'feature':feature,
            'evid':evid,
            'remark':remark,
        })

    logger.info(f"{seqname}: {len(invalid_ratt_features)} RATT features failed validation.")
    ratt_contig_features = get_ordered_features(valid_features)

    logger.info(f"{seqname}: Checking for gene fusion signatures in RATT annotations...")
    ratt_contig_features, merged_features, inconsistent_ratt_features = fusionfisher(
        ratt_contig_features,
        ref_annotation,
    )
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
    evid = None
    remark = None

    if enforce_thresholds:
        ratt_seq_ident = seq_ident
        ratt_seq_covg = seq_covg
    else:
        ratt_seq_ident = ratt_seq_covg = 0

    if feature.location is None:
        valid = False
        evid = "no_coordinates"
        if 'locus_tag' not in feature.qualifiers:
            feature.qualifiers['locus_tag'] = feature.id
        return valid, feature, evid, remark

    if feature.type != 'CDS':
        if feature.type in [
                'rRNA', # these aren't reliably transferred by RATT
                'tRNA',
        ]:
            valid = False
            evid = "categorical"
            remark = f"{feature.type}s categorically rejected in favor of ab initio"
        else:
            valid = True
        return valid, feature, evid, remark

    compound_interval = isinstance(feature.location,Bio.SeqFeature.CompoundLocation)
    # Identify features with 'joins'
    if compound_interval and 'ribosomal_slippage' not in feature.qualifiers:
        #Need to initialize the feature without the compound location attribute.
        #The earliest start and the latest end of the joined feature will be bridged together
        feature_start = feature.location.start
        feature_end = feature.location.end
        feature_strand = feature.location.strand
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
                fix_stop=True,
            )
            if not good_stop:
                evid = "misplaced"
                remark = "RATT-introduced compound interval did not include reference stop position."
                valid = False
                return valid, feature, evid, remark
            # elif feature.corr_possible:
            #     # TODO: find a way to report the stop corrections from the RATT joins without interfering with pseudoscan

            # The stop-corrected coordinates are our new original. pseudoscan will take care of start correction
            # and determining whether that is appropriate to keep.
            feature.og = FeatureProperties()
            feature.corr = FeatureProperties()
            feature.corr_possible = None


    # RATT CDSs don't come with any translations
    feature.qualifiers['translation'] =[str(translate(
        feature.extract(),
        table=cnf.genetic_code,
        to_stop=True,
    ))]

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
        feature_sequence = feature.qualifiers['translation'][0]
        # TODO: if we want to keep track of the blast stats, we should add it as an attribute to the AutarkicSeqFeature object
        #ratt_blast_results.update(blast_stats)

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
                    table=cnf.genetic_code, to_stop=True
                )

            if len(feature_sequence) == 0:
                evid = "zero_length"
                remark = 'length of AA sequence is 0'
            else:
                top_hit, low_covg, blast_stats = BLAST.reference_match(
                    query=SeqRecord(Seq(feature_sequence)),
                    subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                    seq_ident=seq_ident,
                    seq_covg=seq_covg,
                )

                if top_hit:
                    valid = True
                else:
                    evid = "poor_match"
                    remark = 'No blastp hit to corresponding reference CDS at specified thresholds.'


    return valid, feature, evid, remark
