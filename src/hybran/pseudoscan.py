from copy import deepcopy

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import (
    BLAST,
    config,
    extractor,
    designator,
)
from .config import cnf
from .demarcate import (
    coord_check,
    has_broken_stop,
    has_valid_start,
)


def divisible_by_three(feature):
    return len(feature.location) % 3 == 0

def reset_pseudo(feature):
    """
    Remove pseudoscan notes and the pseudo attribute. Useful if reprocessing.
    :param feature: SeqFeature to reset
    :returns feature: reference to the input SeqFeature (which is modified directly)
    """
    keeper_notes = []
    if 'note' in feature.qualifiers:
        keeper_notes = [_ for _ in feature.qualifiers['note'] if not _.startswith('Hybran/Pseudoscan')]
        feature.qualifiers['note'] = keeper_notes

    feature.ps_evid = []
    feature.qualifiers.pop('pseudo', None)
    feature.qualifiers.pop('pseudogene', None)
    return feature

def pseudoscan(
        feature,
        ref_feature,
        seq_ident,
        seq_covg,
        attempt_rescue=False,
        blast_hit_dict=None,
):
    """
    Determine whether a feature should have the /pseudo qualifier.

    :param feature: SeqFeature object of the one to test for pseudo
    :param ref_feature: SeqFeature object of the reference feature to compare to
    :param seq_ident:
    :param seq_covg:
    :param attempt_rescue: Boolean whether to attempt coordinate correction (feature may still be pseudo after correction)
    :param blast_hit_dict:
        dictionary of blast scores for the reference comparison.
        This will be updated if the coordinates are corrected.
    """
    # Check if the reference annotation was pseudo (some branches of the upcoming conditional use this information)
    # seqret moves the reference pseudo tags to a note field beginning with *pseudo or *pseudogene
    # Save the information and drop the note.

    pseudo_note = False
    if 'note' in feature.qualifiers:
        pseudo_note = [_ for _ in feature.qualifiers['note'] if _.startswith("*pseudo") or "Reference gene is pseudo" in _]
        if pseudo_note:
            feature.qualifiers['note'].remove(pseudo_note[0])
        # We're about to decide for ourselves whether the gene is pseudo in this current run
        # and we don't want the existence of the qualifier from the previous run
        # confounding the ref_was_pseudo determination.
        reset_pseudo(feature)

    # checking the feature's pseudo attribute to decide whether the reference is pseudo
    # is only possible with RATT annotations since liftover from the reference may bring
    # the reference's pseudo tag attribute along with it. Our liftover doesn't do that.
    ref_was_pseudo = pseudo_note or designator.is_pseudo(ref_feature.qualifiers)
    ref_d3 = divisible_by_three(ref_feature)
    ref_len = len(ref_feature.location)
    og_feature = deepcopy(feature)

    (feature.og.rcs, feature.og.rce) = coord_check(feature, ref_feature)
    coords_ok = [feature.og.rcs, feature.og.rce]
    og_broken_stop, og_stop_note = has_broken_stop(feature)
    feature.og.d3 = divisible_by_three(feature)
    feature.og.vs = has_valid_start(feature)
    feature.og.ve = not og_broken_stop

    if ref_was_pseudo:
        new_note = []

        broken_stop, stop_note = og_broken_stop, og_stop_note
        if "No stop codons detected" in stop_note:
            (feature.corr.rcs, feature.corr.rce) = coord_check(feature, ref_feature, seek_stop=True)
            feature.corr_accepted = True
            broken_stop, stop_note = has_broken_stop(feature)
            feature.corr.d3 = divisible_by_three(feature)
            feature.corr.vs = has_valid_start(feature)
            feature.corr.ve = not broken_stop

    else:
        fix_start = False
        fix_stop = False
        confirmed_feature = False
        broken_stop, stop_note = has_broken_stop(feature)

        while not confirmed_feature:
            if fix_start or fix_stop:
                (feature.corr.rcs, feature.corr.rce) = coord_check(feature, ref_feature, fix_start, fix_stop)
                coords_ok = [feature.corr.rcs, feature.corr.rce]
                broken_stop, stop_note = has_broken_stop(feature)

                feature.corr.d3 = divisible_by_three(feature)
                feature.corr.vs = has_valid_start(feature)
                feature.corr.ve = not broken_stop

            ref_seq = ref_feature.qualifiers['translation'][0]
            feature_seq = feature.qualifiers['translation'][0]

            #ref_match with 'thresholds enforced'
            top_hit, low_covg, blast_stats = BLAST.reference_match(
                query=SeqRecord(Seq(feature_seq)),
                subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                seq_ident=seq_ident,
                seq_covg=seq_covg,
            )
            blast_ok = top_hit and not low_covg

            if feature.og.bok is None:
                feature.og.bok = bool(blast_ok)
                og_blast_stats = blast_stats
                if all([feature.og.rcs, feature.og.rce]):
                    confirmed_feature = True
                elif attempt_rescue:
                    fix_start = True
                    # If a gene has lost its stop codon and extended (gene longer than reference),
                    # it will have a "bad_stop" but we will want to keep that bad stop in order to
                    # capture as much information as we can (how far along the gene extended).
                    # However, if it has an early stop (gene is shorter than reference), we want to find out
                    # where it was "supposed" to stop (according to reference) and set the stop position there
                    # to capture as much information as we can.
                    if len(feature.location) < len(ref_feature.location):
                        fix_stop = True
                    else:
                        fix_stop = False
            else:
                feature.corr.bok = bool(blast_ok)
                lost_match = (
                    (feature.og.bok and not feature.corr.bok) or
                    (has_valid_start(og_feature) and not has_valid_start(feature) and has_valid_start(ref_feature))
                )
                # Revert coordinate correction if the blastp hit is lost OR the feature does not have a valid start codon
                # (usually due to the reference-corresponding start being affected by an early frameshift)
                if lost_match:
                    feature.corr_accepted = False
                    feature.location = og_feature.location
                    feature.qualifiers = og_feature.qualifiers
                    blast_stats = og_blast_stats
                    stop_note = og_stop_note
                elif feature.corr_possible:
                    feature.corr_accepted = True
                confirmed_feature = True

        # feature is now confirmed
        if blast_hit_dict:
            blast_hit_dict.update(blast_stats[ref_feature.qualifiers['gene'][0]])

    #Remove all transl_except qualifiers from the feature here because they are
    #referencing the location in the reference genome.
    if 'transl_except' in feature.qualifiers:
        feature.qualifiers.pop('transl_except', None)

    #If we managed to confirm the existence of a valid selenocysteine in coord_check,
    #the transl_except attribute will have stored the updated transl_except qualifier string/location.
    if feature.transl_except:
        designator.append_qualifier(
            feature.qualifiers,
            'transl_except',
            feature.transl_except,
        )
    return call(feature, ref_was_pseudo, ref_d3, ref_len, seq_ident, seq_covg)

def call(
        feature,
        ref_was_pseudo,
        ref_d3=None,
        ref_len=None,
        blast_seq_ident=None,
        blast_seq_covg=None,
):
    """
    Determine whether feature is pseudo based on its FeatureProperties.
    Modify the input feature to apply tags and notes accordingly.

    :param feature: AutarkicSeqFeature with populated FeatureProperties
    :param ref_was_pseudo: bool
        Whether the feature's reference was pseudo, as this affects the criteria
    :param ref_d3: bool
        Whether the reference feature is divisible by 3.
        Only applicable if ref_was_pseudo == True.
    :param ref_len: int length of reference feature (for noting length differences)
    :param blast_seq_ident: float alignment sequence identity threshold used
    :param blast_seq_covg: float alignment coverage threshold used
    :returns: True if pseudo / False otherwise.
    """

    #Assign feature properties
    #points to '.corr' if corr_accepted == True
    #points to '.og. if corr_accepted == False
    coords_ok = [feature.rcs, feature.rce]
    # broken_stop is stored in feature.ve, but the note isn't
    broken_stop, stop_note = has_broken_stop(feature)
    d3 = feature.d3 = divisible_by_three(feature)
    valid_start = feature.vs = has_valid_start(feature)
    feature.ve = not broken_stop
    blast_ok = feature.bok

    have_blast_info = (
        blast_ok is not None
        and blast_seq_ident is not None
        and blast_seq_covg is not None
    )

    # Notes that are only interesting if we end up tagging the gene a certain way.
    new_note = []

    if ref_was_pseudo:
        coord_note = (
            f"Locus {'has' if all(coords_ok) else 'does not have'} reference-corresponding "
            f"{'start' if not coords_ok[0] else ''}"
            f"{' and ' if not any(coords_ok) else ''}"
            f"{'end' if not coords_ok[1] else ''}"
            f"{'start and end' if all(coords_ok) else ''}"
        )
        feat_div_note = (
            f"Locus has {'valid' if d3 else 'invalid'} reading frame"
            f"{'' if feature.d3 else '-- not divisible by three'}"
        )
        ref_div_note = (
            f"Reference gene has {'valid' if ref_d3 else 'invalid'} reading frame"
            f"{'' if ref_d3 else '-- not divisible by three'}"
        )
        broke_note = f"{'No internal stop codons and ends with a valid stop codon' if not broken_stop else stop_note}"
        if feature.de:
            new_note.append("Locus has a delayed stop codon")

        if (not all(coords_ok)) and (d3 and not ref_d3) and not broken_stop:
            feature.qualifiers.pop('pseudo', None)
            feature.qualifiers.pop('pseudogene', None)
            is_pseudo = False
        else:
            feature.qualifiers['pseudo'] = ['']
            is_pseudo = True
        feature.ps_evid.append("ref_pseudo")
        new_note.append(f"Reference gene is pseudo")
        new_note.extend([ref_div_note, feat_div_note, coord_note, broke_note])
        new_note = ' | '.join(new_note)
        designator.append_qualifier(feature.qualifiers, 'note', f'Hybran/Pseudoscan: {new_note}')
    else:
        feature.alts = (not feature.rcs and feature.vs)
        feature.alte = (not feature.rce and feature.ve)

        # Summarize coord_check status
        coord_note = (
            f"Locus {'has' if all(coords_ok) else 'does not have'} reference-corresponding "
            f"{'start' if not coords_ok[0] else ''}"
            f"{' and ' if not any(coords_ok) else ''}"
            f"{'end' if not coords_ok[1] else ''}"
            f"{'start and end' if all(coords_ok) else ''}"
        )
        div_note = (
            f"Locus has {'valid' if d3 else 'invalid'} reading frame"
            f"{'' if d3 else '-- not divisible by three'}"
        )
        if have_blast_info:
            blast_note = (
                f"{'Strong' if blast_ok else 'Poor'} blastp match at "
                f"{blast_seq_ident}% identity and {blast_seq_covg}% coverage thresholds"
            )
        else:
            blast_note = ''
        start_note = f"Locus has {'valid' if valid_start else 'invalid'} start codon"
        broke_note = f"{'No internal stop codons and ends with a valid stop codon' if not broken_stop else stop_note}"

        #Note codes are categorized by a '0' or '1' and correspond to 'False' and 'True' respectively
        #D3 = Divisible by three [0/1]
        #VS = Valid start [0/1]
        #VE = Valid end [0/1]
        #RCS = Reference corresponding start [0/1]
        #RCE = Reference corresponding end [0/1]
        #BOK = Blast OK [0/1]
        # represent boolean values as ints but account for N/A (None)
        nacast = lambda _: int(_) if _ is not None else "."
        note_codes = ';'.join([
            f"D3{nacast(feature.d3)}",
            f"VS{nacast(feature.vs)}",
            f"VE{nacast(feature.ve)}",
            f"RCS{nacast(feature.rcs)}",
            f"RCE{nacast(feature.rce)}",
            f"BOK{nacast(feature.bok)}",
        ])
        #This is the code for a normal non-pseudo gene, with no interesting characteristics.
        #These codes will not be added to the feature notes.
        if note_codes == f'D31;VS1;VE1;RCS1;RCE1;BOK{"1" if have_blast_info else "."}':
            note_codes = None

        is_pseudo = True
        if d3 and not broken_stop:
            if all(coords_ok):
                is_pseudo = False
            elif blast_ok:
                is_pseudo = False
            else:
                is_pseudo = True

        if not is_pseudo:
            #Uninteresting cases that do not warrant a 'pseudo note'. Everything is good.
            if all(coords_ok) and blast_ok:
                return is_pseudo

            #Either all(coords_ok) or blast_ok will be True, but not both.
            else:
                #Primary reason for non-pseudo is strong blastp. Interesting because all(coord_ok) == False.
                if blast_ok:
                    new_note.extend([blast_note, coord_note])
                    new_note.extend([broke_note, div_note])

                    if feature.alts:
                        new_note.append("Locus has a valid alternative start site")
                        feature.ps_evid.append("alt_start")

                    if feature.alte:
                        new_note.append("Locus has a valid alternative stop codon")
                        feature.ps_evid.append("alt_end")

                    if feature.de:
                        new_note.append("Locus has a delayed stop codon")

                #Primary reason for non-pseudo is all(coords_ok). Interesting because blast_ok == False.
                elif blast_ok is not None:
                    new_note.extend([coord_note, blast_note])
                    new_note.extend([broke_note, div_note])
                    feature.ps_evid.append('noisy_seq')

                    #Can only comment on differences in gene length if all(coords_ok).
                    if ref_len is not None:
                        if ref_len > len(feature):
                            new_note.append(f"Locus is {ref_len - len(feature)} base pair(s) shorter than the reference")
                        elif ref_len < len(feature):
                            new_note.append(f"Locus is {len(feature) - ref_len} base pair(s) longer than the reference")

        #is_pseudo == True
        #Either not d3 or has broken_stop OR Both all(coords_ok) and blast_ok == False
        else:
            feature.qualifiers['pseudo']=['']

            #Primary reason for pseudo is both all(coords_ok) and blast_ok == False.
            if d3 and not broken_stop:
                new_note.extend([coord_note, blast_note])
                new_note.extend([broke_note, div_note])
                feature.ps_evid.append('no_rcc')
            #Primary reason for pseudo is broken stop or invalid reading frame.
            else:
                if not d3:
                    feature.ps_evid.append('not_div_by_3')
                if broken_stop:
                    feature.ps_evid.append('internal_stop')

                new_note.extend([broke_note, div_note])
                new_note.extend([coord_note, blast_note])

            if not all(coords_ok):
                if feature.alts:
                    new_note.append("Locus has a valid alternative start site")

                if feature.alte:
                    new_note.append("Locus has a valid alternative stop codon")

                if feature.de:
                    new_note.append("Locus has a delayed stop codon")
            else:
                #Can only comment on differences in gene length is all(coords_ok).
                if ref_len is not None:
                    if ref_len > len(feature):
                        new_note.append(f"Locus is {ref_len - len(feature)} base pair(s) shorter than the reference")
                    elif ref_len < len(feature):
                        new_note.append(f"Locus is {len(feature) - ref_len} base pair(s) longer than the reference")

        if new_note:
            new_note = ' | '.join(filter(None, new_note))
            designator.append_qualifier(
                feature.qualifiers,
                'note',
                f'Hybran/Pseudoscan:description:{new_note}',
            )

        if feature.ps_evid:
            designator.append_qualifier(
                feature.qualifiers,
                'note',
                f'Hybran/Pseudoscan:evidence:{";".join(feature.ps_evid)}',
            )

        if note_codes:
            designator.append_qualifier(
                feature.qualifiers,
                'note',
                f'Hybran/Pseudoscan:barcode:{note_codes}',
            )

    return is_pseudo
