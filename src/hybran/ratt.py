import logging

from Bio import SeqIO

from .annomerge import get_and_remove_ref_tracer
from .annomerge import fusionfisher
from . import designator


def postprocess(
        ratt_gbk,
        seq_ident,
        seq_covg,
        nproc=1,
        enforce_thresholds=False,
):
    ratt_contig_record = SeqIO.read(ratt_gbk, 'genbank')
    global record_sequence
    record_sequence = ratt_contig_record.seq

    ratt_contig_non_cds = []
    for feature in ratt_contig_features:
        ref_contig_id = get_and_remove_ref_tracer(feature)
        feature.source = ref_contig_id
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

    logger.debug(f'{seqname}: {len(ratt_contig_non_cds)} non-CDS elements')

    ratt_contig_features, ratt_blast_results, invalid_ratt_features = \
    isolate_valid_ratt_annotations(feature_list=ratt_contig_features,
                                   seq_ident=seq_ident,
                                   seq_covg=seq_covg,
                                   ratt_enforce_thresholds=ratt_enforce_thresholds,
                                   nproc=nproc,
    )
    logger.info(f"{seqname}: {len(invalid_ratt_features)} RATT features failed validation.")
    ratt_contig_features = get_ordered_features(ratt_contig_features)

    logger.info(f"{seqname}: Checking for gene fusion signatures in RATT annotations...")
    ratt_contig_features, merged_features, inconsistent_ratt_features = fusionfisher(ratt_contig_features)
    logger.info(f"{seqname}: {len(merged_features)} RATT fusion genes detected.")
    logger.info(f"{seqname}: {len(inconsistent_ratt_features)} RATT annotations found to be inconsistent.")

    ratt_rejects += invalid_ratt_features + inconsistent_ratt_features

    if merged_features:
        merged_features_record = ratt_contig_record[:]
        merged_features_record.features = merged_features
        SeqIO.write(merged_features_record, output_merged_genes, 'genbank')

    return ratt_contig_features

def isolate_valid_ratt_annotations(feature_list, seq_ident, seq_covg, ratt_enforce_thresholds,
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

    def refcheck(cds_feature, ratt_seq_ident=ratt_seq_ident, ratt_seq_covg=ratt_seq_covg, record_sequence=record_sequence):
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
        feature_sequence = translate(cds_feature.extract(record_sequence), table=genetic_code, to_stop=True)
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
            feature.location = (FeatureLocation(feature_start, feature_end, strand=feature_strand))
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
