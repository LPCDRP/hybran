import os
import logging

import Bio
from Bio import SeqIO
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
from .annomerge import log_feature_fate
from .annomerge import pseudoscan
from .bio import AutarkicSeqFeature


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

        correction_files = [cf for cf in os.listdir(ratt_outdir) if cf.endswith('.Report.txt')]
        for corr_file in correction_files:
            corr_file_path = ratt_file_path + '/' + corr_file
            ratt_correction_files.append(corr_file_path)
    except OSError:
        logger.error('Expecting RATT annotation files but found none')
    if not ratt_features:
        logger.error('RATT did not complete running. Please see the log for more details.')


    with open(invalid_features_logfile, 'w') as rejects_log:
        [log_feature_fate(_[0], rejects_log, _[1]) for _ in invalid_features]

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
    ratt_contig_non_cds = []
    for feature in ratt_contig_record.features:
        feature = AutarkicSeqFeature.fromSeqFeature(feature)
        ref_contig_id = get_and_remove_ref_tracer(feature)
        feature.source = ref_contig_id
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
                'CDS',
                'gene',
                'rRNA', # these aren't reliably transferred by RATT
                'tRNA',
        ]:
            ratt_contig_non_cds.append(feature)
        else:
            ratt_contig_features.append(feature)


    logger.debug(f'{seqname}: {len(ratt_contig_non_cds)} non-CDS elements')

    ratt_contig_features, ratt_blast_results, invalid_ratt_features = \
    isolate_valid_ratt_annotations(
        feature_list=ratt_contig_features,
        ref_annotation=ref_annotation,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
        ratt_enforce_thresholds=enforce_thresholds,
        nproc=nproc,
    )
    logger.info(f"{seqname}: {len(invalid_ratt_features)} RATT features failed validation.")
    ratt_contig_features = get_ordered_features(ratt_contig_features)

    logger.info(f"{seqname}: Checking for gene fusion signatures in RATT annotations...")
    ratt_contig_features, merged_features, inconsistent_ratt_features = fusionfisher(ratt_contig_features)
    logger.info(f"{seqname}: {len(merged_features)} RATT fusion genes detected.")
    logger.info(f"{seqname}: {len(inconsistent_ratt_features)} RATT annotations found to be inconsistent.")

    ratt_rejects = invalid_ratt_features + inconsistent_ratt_features
    ratt_contig_record.features = ratt_contig_features

    # if merged_features:
    #     merged_features_record = ratt_contig_record[:]
    #     merged_features_record.features = merged_features
    #     SeqIO.write(merged_features_record, output_merged_genes, 'genbank')

    return ratt_contig_record, ratt_rejects

def isolate_valid_ratt_annotations(
        feature_list,
        ref_annotation,
        seq_ident,
        seq_covg,
        ratt_enforce_thresholds,
        nproc=1,
):
    """
    This function takes as input a list of features and checks if the length of the CDSs are divisible by
    3 and if the CDS is split across multiple locations. If so, it outputs the features to stdout and removes them
    from the valid_ratt_annotations. The function BLASTs the sequence to corresponding amino acid sequence in
    reference as well.
    :param feature_list: List of features from RATT (list of SeqFeature objects)
    :return: List of valid RATT features (list of SeqFeature objects)
    """
    logger = logging.getLogger('ValidateRATTCDSs')
    logger.debug('Parsing through RATT annotations')
    unbroken_cds = []
    non_cds_features = []
    ratt_blast_results = {}
    rejects = []
    valid_features = []

    if ratt_enforce_thresholds:
        ratt_seq_ident = seq_ident
        ratt_seq_covg = seq_covg
    else:
        ratt_seq_ident = ratt_seq_covg = 0

    def refcheck(cds_feature, ratt_seq_ident=ratt_seq_ident, ratt_seq_covg=ratt_seq_covg):
        valid = False
        remark = ''
        blast_stats = {}
        ref_feature = ref_annotation[
            key_ref_gene(cds_feature.source, cds_feature.qualifiers['gene'][0])
        ]
        try:
            ref_seq = ref_feature.qualifiers['translation'][0]
        except KeyError:
            ref_seq = translate(
                extractor.get_seq(ref_feature),
                table=genetic_code, to_stop=True
            )
        feature_sequence = translate(cds_feature.extract(), table=genetic_code, to_stop=True)
        if len(feature_sequence) == 0:
            remark = 'length of AA sequence is 0'
        else:
            top_hit, low_covg, blast_stats = BLAST.reference_match(
                query=SeqRecord(feature_sequence),
                subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                seq_ident=ratt_seq_ident,
                seq_covg=ratt_seq_covg,
            )

            if top_hit:
                valid = True
            else:
                remark = 'No blastp hit to corresponding reference CDS at specified thresholds.'
        return valid, feature_sequence, blast_stats, remark

    for feature in feature_list:
        if feature.type != 'CDS':
            non_cds_features.append(feature)
            continue

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
                    rejects.append((feature, "RATT-introduced compound interval did not include reference " +
                                    "stop position."))
                    continue

        feature_is_pseudo = pseudoscan(
            feature,
            ref_annotation[key_ref_gene(feature.source, feature.qualifiers['gene'][0])],
            seq_ident,
            seq_covg,
            attempt_rescue=True
        )

        if feature_is_pseudo:
            valid_features.append(feature)
        else:
            unbroken_cds.append(feature)

    logger.debug("Valid CDSs before checking coverage: " + str(len(unbroken_cds) + len(valid_features)))
    logger.debug(f"Checking similarity to reference CDSs using {nproc} process(es)")

    with multiprocessing.Pool(processes=nproc) as pool:
         results = pool.map(
            refcheck,
            unbroken_cds,
        )
    for i in range(len(unbroken_cds)):
        valid, feature_sequence, blast_stats, rejection_note = results[i]
        cds_feature = unbroken_cds[i]
        if valid:
            ratt_blast_results.update(blast_stats)
            cds_feature.qualifiers['translation'] = [str(feature_sequence)]
            valid_features.append(cds_feature)
        else:
            rejects.append((cds_feature, rejection_note))
    logger.debug("Valid CDSs after checking coverage: " + str(len(valid_features)))
    return valid_features, ratt_blast_results, rejects
