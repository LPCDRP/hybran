import logging
from copy import deepcopy
from math import log, ceil

from Bio import Align
from Bio.Data import CodonTable
from Bio.Seq import translate
from Bio.SeqFeature import (
    FeatureLocation,
    ExactPosition,
    SeqFeature,
    BeforePosition,
    AfterPosition,
)

from . import (
    config,
    designator,
    extractor,
)
from .config import cnf


def has_valid_start(feature):
    """
    Finds the first three base pairs in a SeqFeaure and determines if it
    has a valid start codon according to its corresponding codon table.
    :param feature: A SeqFeature object
    :return: True if valid start codon exists
    """
    start_codons = CodonTable.generic_by_id[cnf.genetic_code].start_codons
    feature_seq = str(feature.extract(parent_sequence=feature.references[feature.location.parts[0].ref],
                                      references=feature.references))[:3]

    return feature_seq in start_codons

def has_broken_stop(feature):
    """
    Finds the amount and location of internal stops.
    :param feature: A SeqFeature object
    """
    internal_stop = False
    note = ''
    feature_seq = feature.extract(parent_sequence=feature.references[feature.location.parts[0].ref],
                                  references=feature.references)
    translation = str(feature_seq.translate(to_stop=False, table=cnf.genetic_code))
    num_stop = [i for i,e in enumerate(translation) if e == "*"]
    num_internal_stop = [i for i,e in enumerate(translation) if e == "*" and i != (len(translation)-1)]
    if len(num_internal_stop) >= 1 or translation[-1] != "*":
        internal_stop = True
        if len(num_stop) == 0:
            note = f"No stop codons detected in the translated sequence"
        else:
            note = f"Internal stop detected at codon(s) {' '.join([str(i) for i in num_internal_stop])}"
    return internal_stop, note

def stopseeker(feature, circularize=False):
    """
    Use the coordinates from a feature that does not contain a valid in-frame stop codon, and
    return a new extended feature that contains the next available downstream stop codon.
    :param feature: A SeqFeature object
    :param circular: User specified Boolean to determine topology (i.e. 'circular' or 'linear')
    :return: SeqFeature object that contains a valid stop codon
    """
    extended_feature = deepcopy(feature)
    feature_ref = feature.location.parts[0].ref
    record_sequence = feature.references[feature.location.parts[0].ref]
    rec_len = len(record_sequence) -1

    if feature.strand == 1:
        extended_feature.location.parts[-1]._end = ExactPosition(len(record_sequence) - 1)
    else:
        extended_feature.location.parts[0]._start = ExactPosition(0)

    if circularize and ((0 in feature) or (rec_len in feature)):
        circ = FeatureLocation(int(0), int(rec_len), strand=feature.strand, ref=feature_ref)
        #Creates CompoundLocation
        if feature.strand == 1:
            extended_feature.location = extended_feature.location + circularize
        else:
            extended_feature.location = circularize + extended_feature.location
    extended_seq = extended_feature.extract().translate(to_stop=True, table=cnf.genetic_code)
    #Translation extends up to, but not including the stop codon, so we need to add 3 to the end
    extended_seq_len = len(extended_seq)*3 + 3
    original_seq_len = len(feature.extract())
    len_difference = extended_seq_len - original_seq_len

    if feature.strand == 1:
        extended_feature_end = feature.location.parts[-1].end + len_difference
        if extended_feature_end > rec_len:
            #Annotation will look like [feature.start ... >(rec_len)]
            extended_feature_end = AfterPosition(rec_len)

        extended_feature.location.parts[-1]._end = ExactPosition(
            extended_feature_end
        )
    else:
        extended_feature_start = feature.location.parts[0].start - (len_difference - 1)
        if extended_feature_start < 0:
            #Annotation will look like [<0 ... feature.end]
            extended_feature_start = BeforePosition(0)

        extended_feature.location.parts[0]._start = ExactPosition(
            extended_feature_start
        )

    feature.location = extended_feature.location

    return feature

def update_termini(base_location, start, end):
    """
    Update a feature location object with a new start and end position,
    while preserving compound location components if the exist.
    :param base_location: FeatureLocation or CompoundLocation object
    :param start: int or ExactPosition or the like new start position
    :param end: int or ExactPosition or the like new end position
    :returns: modified input location object. the given object is also modified directly
    """
    #
    # Setting .start and .end directly is not allowed, but these
    # _start and _end properties have the right effect.
    # If it ever breaks, the recourse will be to make a new location object
    # (either SimpleLocation or CompoundLocation) and replicate everything
    # but the parts[0] start and parts[-1] end.
    #
    base_location.parts[0]._start = ExactPosition(start)
    base_location.parts[-1]._end = ExactPosition(end)
    return base_location

def coord_check(
        feature,
        ref_feature,
        fix_start=False,
        fix_stop=False,
        seek_stop=None,
        ref_gene_name=None,
):
    """
    This function takes a feature as an input and aligns it to the corresponding reference gene.
    This function will return two Boolean values (True if the start/end coordinates align with the reference).
    The start/end coordinates can be corrected to match the reference if desired in various circumstances.
    :param feature: SeqFeature object
    :param fix_start: Boolean
    :param fix_stop: Boolean
    :param seek_stop: Boolean whether to look for a valid stop codon if post-correction.
    :param ref_gene_name: str reference gene to check against.
        Must be a key existing in the `ref_annotation` dictionary.
        If not defined, the reference gene matching `feature`'s gene qualifier is used instead.
    :return: True/False if the start/stop was fixed
    """

    logger = logging.getLogger('CoordCheck')
    record_sequence = feature.references[feature.location.parts[0].ref]
    ref_seq = extractor.get_seq(ref_feature)
    feature_start = feature.location.start
    feature_end = feature.location.end
    feature_seq = feature.extract()
    og_feature = deepcopy(feature)
    if 'gene' not in og_feature.qualifiers:
        og_feature.qualifiers['gene'] = [ref_gene_name]
    og_feature_start = og_feature.location.start
    og_feature_end = og_feature.location.end
    og_feature_loc_ref = og_feature.location.parts[0].ref

    if feature.og.location is None:
        feature.og.location = og_feature.location

    def coord_align(ref_seq, feature_seq):
        """
        This function generates an alignment between a feature and a reference sequence.
        :param ref_seq:
        :param feature_seq:
        :returns:
            -found_low - Boolean - if the feature sequence contains a reference corresponding start
            -found_high - Boolean - if the feature sequence contains a reference corresponding end
            -target - A numpy.ndarray of aligned positions with respect to reference sequence
            -query - A numpy.ndarray of aligned positions with respect to the feature sequence
            -alignment - A Bio.Align.Alignment object illustrating the pairwise sequence alignment
            -padding - Boolean - if the feature sequence is capable of adding up/downstream context (padding)
            -score - The alignment score (float)
            -interval - The number of consecutive matching base pairs needed to establish reference correspondence
            -relaxed_found_high - Boolean - feature contains a reference corresponding stop (ignoring mismatches)
        """
        #Probability that the continuous interval used to find good start/stops
        #occurs by chance should be = 1/(ref_seq*10)
        interval = max(ceil((log(len(ref_seq) * 10))/log(4)), 3)

        aligner = Align.PairwiseAligner(scoring="blastn", mode = 'global')
        #"blastn" scoring: extend_gap_score = -2.0, open_gap_score = -7.0
        #Punishing internal gaps slighly more will discourage aberrant behavior from the aligner.
        #Prevents the miscalling of found_low/high is some scenarios and encourages continuity.
        #Ex)
        #internal_open_gap_score= -7.0           internal_open_gap_score= -8.0
        #internal_extend_gap_score = -2.0        internal_extend_gap_score = -2.01
        # GT--------AG                           GTAG--------
        # ||--------||                   ->      ||||--------
        # GTAGTCCGCGAG                           GTAGTCCGCGAG
        aligner.internal_open_gap_score = -8.0
        aligner.internal_extend_gap_score = -2.01

        alignment = aligner.align(ref_seq, feature_seq)

        alignment = alignment[0]
        target = alignment.aligned[0]
        query = alignment.aligned[1]
        score = alignment.score
        padding = False

        found_low = (target[0][0] == 0) and (abs(target[0][0] - target[0][1])) >= interval
        found_high = (target[-1][1] == len(ref_seq)) and (abs(target[-1][0] - target[-1][1])) >= interval

        #Sequence intervals of the lowest and highest aligned bases.
        target_low_seq = get_gapped_sequence(alignment, 'target', target[0][0], target[0][1])
        target_high_seq = get_gapped_sequence(alignment, 'target', target[-1][0], target[-1][1])
        query_low_seq = get_gapped_sequence(alignment, 'query', query[0][0], query[0][1])
        query_high_seq = get_gapped_sequence(alignment, 'query', query[-1][0], query[-1][1])

        relaxed_found_high = found_high
        if len(target) > 1 or len(query) > 1: #more than one interval blocks exist in the alignment sequence - discontinuous at some point.
            target_inter_gaps = get_gapped_sequence(alignment, 'target', target[-2][0], target[-1][1]).count("-")
            query_inter_gaps = get_gapped_sequence(alignment, 'query', query[-2][0], query[-1][1]).count("-")

            target_penultimate_interval_seq = get_gapped_sequence(alignment, 'target', target[0][0], target[-2][1])[-interval:]
            query_penultimate_interval_seq = get_gapped_sequence(alignment, 'query', query[0][0], query[-2][1])[-interval:]

            penultimate_found_high = (
                # 1) Penultimate interval need to match exactly (last ~7 bp of second to last alignment block)
                (target_penultimate_interval_seq == query_penultimate_interval_seq)
                # 2) The entire penultimate alignment block needs to be greater than the interval (~7)
                and (abs(target[-2][0] - target[-2][1]) >= interval)
                # 3) The number of gaps between the last and second to last alignment blocks cannot exceed the interval (~7)
                #and (target_inter_gaps + query_inter_gaps <= interval)
            )
            relaxed_found_high = penultimate_found_high or found_high

        #make sure there aren't too many mismatches causing falsely assigned found_low/high values
        #
        # relaxed_found_high was determined using found_high before the following adjustments
        # because we need it to be True in the case of non-stop SNPs. The final found_high
        # in these cases should be false, though, so the following achieves that.
        if (
                (target_low_seq[:3] != query_low_seq[:3])
                and not ([has_valid_start(ref_feature), has_valid_start(feature)] == [True, True])
        ):
            found_low = False
        elif found_low and target[0][1] < (len(ref_seq)/3):
            gaps = ident = mismatch = 0
            for i in range(len(target_low_seq)):
                if target_low_seq[i] == '-' or query_low_seq[i] == '-':
                    gaps += 1
                elif target_low_seq[i] == query_low_seq[i]:
                    ident += 1
                else:
                    mismatch += 1
            if (gaps + mismatch) > (len(target_low_seq)/3):
                found_low = False

        if (target_high_seq[-3:] != query_high_seq[-3:]):
            found_high = False
        elif found_high and abs(target[-1][1] - target[-1][0]) < (len(ref_seq)/3):
            gaps = ident = mismatch = 0
            for i in range(len(target_high_seq)):
                if target_high_seq[i] == '-' or query_high_seq[i] == '-':
                    gaps += 1
                elif target_high_seq[i] == query_high_seq[i]:
                    ident += 1
                else:
                    mismatch += 1
            if (gaps + mismatch) > (len(target_high_seq)/3):
                found_high = False

        if not found_low or not found_high:
            padding = True
        return found_low, found_high, target, query, alignment, padding, score, interval, relaxed_found_high

    def get_gapped_sequence(alignment, seq_type, start, stop):
        seq_types = ['target', 'query']
        if seq_type not in seq_types:
            raise ValueError(f"Invalid sequence type. Expected one of: {seq_types}")
        elif seq_type == 'target':
            gapped_seq = list(alignment.indices[0])
            alignment = alignment[0]
        else:
            gapped_seq = list(alignment.indices[1])
            alignment = alignment[1]

        #The index of the stop position is one off from the stop position itself
        stop -= 1

        interval_seq = alignment[gapped_seq.index(start) : gapped_seq.index(stop) + 1]
        return interval_seq

    def add_padding(feature, target, query, interval):
        """
        If we're looking to make corrections, this function will add up/downstream context
        to help the aligner search for reference correspondence.
        :param feature: An AutarkicSeqFeature object
        :param target: A numpy.ndarray of aligned positions with respect to reference sequence
        :param query: A numpy.ndarray of aligned positions with respect to the feature sequence
        :param interval: The number of consecutive matching base pairs needed to establish reference correspondence
        :return pad_feature: The AutarkicSeqFeature object with additional context (new coordinates).
        """
        pad_feature = deepcopy(feature)
        feature_start = feature.location.start
        feature_end = feature.location.end

        #The only time padding matters is when the reference overhangs on the left or right side of
        #the feature. All other scenarios would be where found_low/high = True and extra context is unnecessary.
        left_ref_overhang = int(target[0][0]) > int(query[0][0])
        right_ref_overhang = (int(query[-1][1]) == len(feature)) and (int(target[-1][1]) < len(ref_seq))
        left_pad = right_pad = abs(len(ref_seq) - len(og_feature.extract())) + 2*interval

        if left_ref_overhang:
            left_pad = int(target[0][0] - query[0][0])
            left_pad = left_pad + 2*interval
        if right_ref_overhang:
            right_pad = len(alignment[0]) - int(query[-1][1])
            right_pad = right_pad + 2*interval

        if (fix_start and feature.strand == 1) or (fix_stop and feature.strand == -1):
            if (fix_stop and feature.strand == -1):
                left_pad = right_pad
            padded_feature_start = max(0, (feature_start - left_pad))
            if (not found_low and feature.strand == 1) or (not found_high and feature.strand == -1):
                feature_start = padded_feature_start

        if (fix_stop and feature.strand == 1) or (fix_start and feature.strand == -1):
            if (fix_start and feature.strand == -1):
                right_pad = left_pad
            padded_feature_end = min(len(record_sequence), feature_end + right_pad)
            if (not found_low and feature.strand == -1) or (not found_high and feature.strand == 1):
                feature_end = padded_feature_end

        update_termini(pad_feature.location, feature_start, feature_end)
        return pad_feature

    def has_delayed_end(feature_seq, query, found_high, relaxed_found_high, rce):
        """
        An AutarkicSeqFeature object with a "delayed end" is determined by coord_check after an alignment
        has been generated. The "delayed end" criteria is used to identify non-stop mutations and
        potential gene fusion events.
        :param feature_seq: String of feature sequence
        :param query: A numpy.ndarray of aligned positions with respect to the feature sequence
        :param found_high: Boolean - feature contains a reference corresponding stop
        :param relaxed_found_high: Boolean - feature contains a reference corresponding stop (ignoring mismatches)
        :param rce: Boolean - feature has a positionally identical reference corresponding end
        :return delayed_end: Boolean  - whether the feature has a delayed end / non-stop mutation.
        """

        delayed_end = (
            not rce
            and (
                found_high or relaxed_found_high
                )
            and (
                # the end of the feature extends beyond the last reference base
                query[-1][1] < len(feature_seq)
                # The implementation of 'relaxed_found_high' and the change to our internal gap score
                # should prevent spurious alignment blocks induced by delayed stops. Negating the need for
                # the statement below (which was used to identify these spurious alignment blocks that did not
                # extend beyond the last reference base.)
                # or (final_target[-1][1] - final_target[-1][0]) < final_interval
                )
            )
        return delayed_end

    #First alignment
    found_low, found_high, target, query, alignment, padding, first_score, interval, relaxed_found_high = coord_align(ref_seq, feature_seq)
    #Assign initial alignment, but don't overwrite it.
    if feature.og.alignment is None:
        feature.og.alignment = alignment

    corrected_feature = deepcopy(feature)
    corrected_feature_start = corrected_feature.location.start
    corrected_feature_end = corrected_feature.location.end

    if padding:
        #Align again after adding padding to the feature sequence if warranted
        pad_feature = add_padding(feature, target, query, interval)
        pad_feature_seq = pad_feature.extract()
        pad_found_low, pad_found_high, pad_target, pad_query, pad_alignment, padding, second_score, second_interval, pad_relaxed_found_high = coord_align(ref_seq, pad_feature_seq)

        #Don't try to fix what isn't broken
        if found_low:
            pad_found_low = False
        if found_high:
            pad_found_high = False
    else:
        second_score = -1
        pad_found_high = False
        pad_found_low = False

    for i in range(2):
        if feature.strand == 1:
            good_stop = found_high
            if pad_found_high and fix_stop:
                corrected_feature_end = (pad_feature.location.start + pad_query[-1][1])
                if (second_score > first_score):
                    feature_end = corrected_feature_end
                    good_stop = True
            elif found_high:
                corrected_feature_end = feature_start + query[-1][1]
                if fix_stop:
                    feature_end = corrected_feature_end

            good_start = found_low
            if pad_found_low and fix_start:
                corrected_feature_start = (pad_feature.location.start + (pad_query[0][0]))
                if (second_score > first_score):
                    feature_start = corrected_feature_start
                    good_start = True
            elif found_low:
                corrected_feature_start = feature_start + query[0][0]
                if fix_start:
                    feature_start = corrected_feature_start

        elif feature.strand == -1:
            good_stop = found_high
            if pad_found_high and fix_stop:
                corrected_feature_start = (pad_feature.location.end - pad_query[-1][1])
                if (second_score > first_score):
                    feature_start = corrected_feature_start
                    good_stop = True
            elif found_high:
                corrected_feature_start = feature_end - query[-1][1]
                if fix_stop:
                    feature_start = corrected_feature_start

            good_start = found_low
            if pad_found_low and fix_start:
                corrected_feature_end = (pad_feature.location.end - pad_query[0][0])
                if (second_score > first_score):
                    feature_end = corrected_feature_end
                    good_start = True
            elif found_low:
                corrected_feature_end = feature_end - query[0][0]
                if fix_start:
                    feature_end = corrected_feature_end

        #Catch corner cases where we try to correct a reference corresponding start PAST the end of the feature.
        #or where we try to correct a reference corresponding stop BEFORE the start of a feature.
        #If the original alignment was really far off, and we found a downstream start/upstream stop
        #simply by chance from the addition of extra padding, we should just leave the gene alone.
        if feature_start > feature_end:
            feature_start = og_feature_start
            feature_end = og_feature_end
            logger.warning(f"Attempted to correct {feature.qualifiers['gene'][0]} with invalid coordinates. Restoring original positions.")
        update_termini(feature.location, feature_start, feature_end)
        feature_seq = feature.extract()

        if i == 1:
            continue

        #Coordinates can fail to get corrected because the second alignment score is worse than the first score
        #only because of excess leftover padding. A third alignment with the corrected coordinates (sans 'padding')
        #ensures we aren't unnecessarily rejecting corrections.
        elif any([pad_found_low, pad_found_high]) and any([fix_start, fix_stop]) and (first_score > second_score):

            #Same corner case 'catch' as the one found above
            if corrected_feature_start > corrected_feature_end:
                 update_termini(
                     feature.location,
                     og_feature_start,
                     og_feature_end,
                 )
                 logger.warning(f"Attempted to correct {feature.qualifiers['gene'][0]} with invalid coordinates. Restoring original positions.")
                 break
            update_termini(
                corrected_feature.location,
                corrected_feature_start,
                corrected_feature_end,
            )
            corrected_feature_seq = corrected_feature.extract()
            cor_low, cor_high, cor_target, cor_query, cor_alignment, cor_padding, third_score, third_interval, cor_relaxed_found_high = coord_align(ref_seq, corrected_feature_seq)

            if (third_score >= first_score):
                second_score = third_score + 1
                continue
            else:
                break
        else:
            break

    #If no stops are detected in the final feature, run stopseeker to find nearest downstream stop.
    if seek_stop or ((fix_start or fix_stop) and seek_stop is None):
        broken_stop, stop_note = has_broken_stop(feature)
        if "No stop codons detected" in stop_note:
            extended_feature = stopseeker(feature)
            if len(extended_feature) > len(feature):
                feature.location = extended_feature.location

    final_feature_seq = feature.extract()
    final_found_low, final_found_high, final_target, final_query, final_alignment, final_padding, final_score, final_interval, final_relaxed_found_high = coord_align(ref_seq, final_feature_seq)


    #Up to this point, good_start/stop would be True if a sequence CONTAINED a reference corresponding start/stop somewhere
    #within its boundaries. This concept was necessary for determining where and when corrections should be made.
    #Now that corrections have been made, we are requiring that good_start/stop should only be True if positionally identical
    #to the start/stop in the reference sequence.
    if good_start:
        if ((feature.strand == 1 and feature.location.start != corrected_feature_start) or
            (feature.strand == -1 and feature.location.end != corrected_feature_end)):
            good_start = False
    if good_stop:
        if ((feature.strand == 1 and feature.location.end != corrected_feature_end) or
            (feature.strand == -1 and feature.location.start != corrected_feature_start)):
            good_stop = False

    if og_feature.location != feature.location:
        feature.qualifiers['translation'] = [
            str(translate(
                feature.extract(),
                table=cnf.genetic_code,
                to_stop=True,
            ))
        ]
        designator.append_qualifier(
            feature.qualifiers, 'inference',
            "COORDINATES:alignment:Hybran"
        )
        #Assign feature.corr attributes if a change was made
        feature.corr.alignment = final_alignment
        feature.corr.location = deepcopy(feature.location)
        feature.corr_possible = True
    elif feature.corr_possible is None and (fix_start or fix_stop or seek_stop):
        feature.corr_possible = False

    if feature.og.de is None:
        feature.og.de = has_delayed_end(og_feature.extract(), query, found_high, relaxed_found_high, good_stop)
    if feature.corr_possible:
        feature.corr.de = has_delayed_end(final_feature_seq, final_query, final_found_high, final_relaxed_found_high, good_stop)

    return good_start, good_stop
