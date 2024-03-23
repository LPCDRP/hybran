# Title: Merge annotation from Prokka for the positions which are not annotated by RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files in case of
# multiple contigs/chromosome annotations and a Genbank file (.gbk) file from Prokka annotation run with a reference and
# a Genbank file (.gbk) file from Prokka annotation run without a reference. The output is an EMBL file with annotation
# predominantly from RATT and the intergenic regions annotated by RATT are filled with Prokka. This script also
# generates a log file to indicate characteristics of the transferred features from Prokka.

import sys
import collections
import os
import tempfile
import pickle
import logging
import time
import re
from copy import deepcopy

# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from . import BLAST
from . import converter
from . import config
from . import designator
from . import extractor
from . import pseudoscan
from . import __version__
from .bio import SeqIO
from .bio import translate
from .bio import FeatureProperties
from .compare import (
    compare,
    have_same_stop,
    overlap_inframe,
)
from .demarcate import (
    coord_check,
    has_broken_stop,
    update_termini,
)
from .lumberjack import log_feature_fates
from .lumberjack import log_coord_corrections
from .lumberjack import log_pseudos
from .util import keydefaultdict, mpbreakpoint

def get_and_remove_ref_tracer(feature):
    """
    Remove the temporary tracer note we added to keep track of where an annotation that RATT placed originated from
    """
    ref_contig_id = ""
    if 'note' in feature.qualifiers:
        marker_note = [_ for _ in feature.qualifiers['note'] if _.startswith("HYBRANSOURCE")][0]
        feature.qualifiers['note'].remove(marker_note)
        if not feature.qualifiers['note']:
            del feature.qualifiers['note']
        ref_contig_id = re.sub(r'\s+', '', marker_note.split(':')[1])

    return ref_contig_id

def ref_fuse(fusion_gene_name):
    """
    Create a dummy SeqFeature to use as a reference for fusion gene coordinate checking/correction.
    Throws a KeyError if one of the constituent names is not in the ref_annotation dictionary or if the fusion gene name is not valid.
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
        ref_feature_i = ref_annotation[key_ref_gene(const_refs[i], const_genes[i])]
        location_sequences.update(ref_feature_i.references)
        location_parts += list(ref_feature_i.location.parts)
    ref_fusion = SeqFeature(
        # Biopython currently doesn't support the 'order' operator for feature.extract()
        Bio.SeqFeature.CompoundLocation(location_parts, operator='join'),
        qualifiers={'locus_tag':fusion_gene_name, 'gene':fusion_gene_name},
    )

    ref_fusion.references = location_sequences
    return ref_fusion

def key_ref_gene(ref_id, gene_name):
    """
    Generate a key for ref_annotation dictionary

    We do this rather than use a nested dictionary because, for gene fusions,
    the reference sequence that each member lifted over from is not necessarily the same.
    """
    return '@@@'.join([ref_id, gene_name])

def loc_index(feature):
    """
    Produce the index tuple we use to track feature conflicts.
    """
    return (
        int(feature.location.start),
        int(feature.location.end),
        int(feature.location.strand),
    )

def get_ordered_features(feature_list):
    """
    This function takes list of features and returns the list sorted by genomic location
    :param feature_list: list of type SeqFeature (Biopython feature) formats
    :return: sorted list of type SeqFeature (Biopython feature) formats, sorted by genomic location
    """

    features_dict = {}
    ordered_features = []
    for feature in feature_list:
        feature_start = int(feature.location.start)
        if feature_start not in features_dict.keys():
            features_dict[feature_start] = [feature]
        else:
            features_dict[feature_start].append(feature)
    ordered_features_dict = collections.OrderedDict(sorted(features_dict.items()))
    for position in ordered_features_dict.keys():
        features_in_position = ordered_features_dict[position]
        for feature in features_in_position:
            ordered_features.append(feature)
    return ordered_features


def generate_feature_dictionary(feature_list):
    """
    This function takes as input a list of features and returns a dictionary with the key as a tuple of
    feature start and stop positions and the value as the feature.
    :param feature_list: List of features (SeqFeature objects)
    :return: sorted dictionary ordered by the genomic position i.e. feature location where key is a tuple
    (feature_start, feature_end, feature strand) and value is the corresponding SeqFeature object
    """
    feature_dict = {}
    for feature in feature_list:
        feature_key = (int(feature.location.start), int(feature.location.end), int(feature.location.strand))
        feature_dict[feature_key] = feature
    sorted_feature_dict = collections.OrderedDict(sorted(feature_dict.items()))
    return sorted_feature_dict

def merge_qualifiers(f1quals, f2quals):
    """
    Combine two qualifier dictionaries, combining lists for qualifiers that may
    have multiple values.
    :param f1quals: dict first SeqFeature's qualifiers
    :param f2quals: dict second SeqFeature's qualifiers
    :return: dict combined qualifier dictionary
    """
    multifields = [
            'note',
            'gene_synonym',
            'experiment',
            'inference',
    ]
    final_qualifiers = deepcopy(f1quals)
    final_qualifiers.update(f2quals)
    for qual in multifields:
        if qual in f1quals.keys() and qual in f2quals.keys():
            final_qualifiers[qual] = list(
                set(f1quals[qual]).union(set(f2quals[qual]))
            )
    return final_qualifiers

def complementary(f1, f2, f1_status=None, f2_status=None):
    """
    Given two (Autarkic)SeqFeatures, check whether their reference-correspondence is complementary.
    That is, one has a good start and bad stop while the other has a bad start and good stop.
    :param f1: SeqFeature
    :param f2: SeqFeature
    :param f1_status: bool 2-tuple (goodstart, goodstop) for f1
    :param f2_status: bool 2-tuple (goodstart, goodstop) for f2
        If None, attempt to use (f1(2).rcs, f1(2).rce). The argument allows flexibility for
        providing coord_check results for f1 or f2 with respect to alternative reference genes.
    :return: bool whether the two features are complementary
    """
    if f1_status is None:
        f1_status = (f1.rcs, f1.rce)
    if f2_status is None:
        f2_status = (f2.rcs, f2.rce)
    return (
        (any(f1_status) and any(f2_status))
        and (int(f1_status[0])+int(f2_status[0]), int(f1_status[1])+int(f2_status[1]))==(1,1)
    )

def fissionfuser(flist, seq_ident, seq_covg):
    """
    Given a list of features ordered by genomic position, identify adjacent gene fragments and combine them into a single feature.
    :param flist: list of SeqFeature objects
    :param seq_ident: sequence identity threshold for BLAST (for pseudo-calling)
    :param seq_covg: alignment coverage threshold for BLAST (for pseudo-calling)
    :returns:
        list of SeqFeature objects to keep (some modified from the original)
        list of annotations that have been merged into their neighbor.
    """
    logger = logging.getLogger('FissionFuser')
    outlist = []
    dropped_ltag_features = []
    last_gene_by_strand = {}
    for feature in flist:
        outlist.append(feature)
        if feature.type != 'CDS':
            continue
        if feature.location.strand in last_gene_by_strand:
            last_gene = last_gene_by_strand[feature.location.strand]
            # we'll be working with the loop variable so we can update the dictionary ahead of time
            last_gene_by_strand[feature.location.strand] = feature
        else:
            last_gene_by_strand[feature.location.strand] = feature
            continue

        last_gene_named = 'gene' in last_gene.qualifiers
        curr_gene_named = 'gene' in feature.qualifiers
        only_one_named = last_gene_named ^ curr_gene_named

        combine = False
        reason = ''

        # Many ab initio CDSs are rejected before this step, so we cannot deduce adjacency on the genome by adjacency on the list.
        # We can use the difference in their locus tag numbers to determine that.
        n_last_gene = int(last_gene.qualifiers['locus_tag'][0].split('_')[1])
        n_curr_gene = int(feature.qualifiers['locus_tag'][0].split('_')[1])
        dist = n_curr_gene - n_last_gene

        #
        # This function runs after coordinate correction, so if two ab initio CDSs now overlap in-frame,
        # it means that one of them had its start position corrected to correspond to that of the other fragment.
        #
        if overlap_inframe(last_gene.location, feature.location):
            combine = True
            reason = 'overlapping_inframe'
        #
        # Don't merge CDSs that are too far apart.
        # We only want to consider those that are directly adjacent, but we have to account for alternating genes on the opposite strand.
        # Using 3 as a cutoff because I do not expect more than one or two overlapping genes/gene fragments on the opposite strand between two fragments of what should be a single gene.
        #
        elif dist > 3:
            continue
        #
        # When neither are named, we have no other way to know whether they should be combined
        #
        elif (not last_gene_named) and (not curr_gene_named):
            continue
        #
        # Check for complementarity
        #
        elif only_one_named or extractor.get_gene(last_gene) == extractor.get_gene(feature):
            if only_one_named and last_gene_named:
                ref_gene_source = last_gene.source
                ref_gene = last_gene.qualifiers['gene'][0]
            else:
                ref_gene_source = feature.source
                ref_gene = feature.qualifiers['gene'][0]

            lg_status = coord_check(
                last_gene,
                ref_annotation[key_ref_gene(ref_gene_source, ref_gene)],
            )
            cg_status = coord_check(
                feature,
                ref_annotation[key_ref_gene(ref_gene_source, ref_gene)],
            )
            if complementary(
                    f1=last_gene,
                    f1_status=lg_status,
                    f2=feature,
                    f2_status=cg_status,
            ):
                combine = True
                reason = 'complementary_fragments'

        #
        # Combine intervals and check validity, aborting otherwise
        #
        if combine:
            new_start = last_gene.location.start
            new_end = feature.location.end
            if only_one_named and last_gene_named:
                new_feature = deepcopy(last_gene)
                dropped_feature = feature
            else:
                new_feature = deepcopy(feature)
                dropped_feature = last_gene
            dropped_feature_name = f"{extractor.get_ltag(dropped_feature)}:{extractor.get_gene(dropped_feature)}"
            new_feature_name = f"{extractor.get_ltag(new_feature)}:{extractor.get_gene(new_feature)}"
            update_termini(new_feature.location, new_start, new_end)
            new_feature.og = FeatureProperties()
            new_feature.corr = FeatureProperties()
            new_feature.corr_accepted = new_feature.corr_possible = None
            new_feature.qualifiers = merge_qualifiers(dropped_feature.qualifiers, new_feature.qualifiers)

            if reason == 'overlapping_inframe':
                confirmed = True
            #
            # Abort combination if one of the genes is non-pseudo. (issue #66)
            #
            elif last_gene_named and not designator.is_pseudo(last_gene.qualifiers):
                confirmed = False
                problem = f"{last_gene.qualifiers['locus_tag'][0]} is not pseudo"
            elif curr_gene_named and not designator.is_pseudo(feature.qualifiers):
                confirmed = False
                problem = f"{feature.qualifiers['locus_tag'][0]} is not pseudo"
            #
            # Last check before pulling the trigger: coordinate verification
            #
            elif not all(coord_check(
                    new_feature,
                    ref_annotation[key_ref_gene(new_feature.source, new_feature.qualifiers['gene'][0])],
                    fix_start=True, fix_stop=True)
            ):
                confirmed = False
                problem = "coordinate verification failed"
            else:
                confirmed = True


            if confirmed:
                dropped_ltag_features.append({
                    'feature':dropped_feature,
                    'evid':reason,
                    'remark':"combined fission fragments",
                    'superior':new_feature,
                })

                # TODO: pseudoscan is the only thing at this time that determines corr_accepted.
                # If it isn't making the correction itself, it won't make that call, so we're
                # temporarily re-doing the correction to make that happen...
                # TODO: This will make it seem like there isn't a correction happening, but it
                # is necessary to prevent pseudoscan from rejecting our initial correction.
                new_feature.og = FeatureProperties()
                new_feature.corr = FeatureProperties()
                new_feature.corr_accepted = new_feature.corr_possible = None
                # Re-call pseudoscan for updated notes
                pseudoscan.pseudoscan(
                    new_feature,
                    ref_annotation[key_ref_gene(new_feature.source, new_feature.qualifiers['gene'][0])],
                    seq_ident=seq_ident,
                    seq_covg=seq_covg,
                    attempt_rescue=True,
                )
                last_gene_by_strand[feature.location.strand] = new_feature
                outlist.remove(last_gene)
                outlist.remove(feature)
                outlist.append(new_feature)
            else:
                logger.debug(f"Did not combine {dropped_feature_name} and {new_feature_name} ({reason}) because {problem}.")

    return outlist, dropped_ltag_features


def fusion_name(feature1, feature2):
    """
    For use by fusionfisher.
    Create a fusion gene and product name. Usually, it's just "gene1::gene2" and "product1 / product2",
    but we need to watch out for redundancies, especially when one of the two features is already a fusion
    including the second's name.
    """
    f1_name_components = extractor.get_gene(feature1).split('::')
    f2_name_components = extractor.get_gene(feature2).split('::')
    reduced_f2_name_components = [_ for _ in f2_name_components if _ not in f1_name_components]
    name = '::'.join(f1_name_components + reduced_f2_name_components)

    f1_source_components = feature1.source.split('::')
    f2_source_components = feature2.source.split('::')
    # only get the host sequence IDs corresponding to the gene names we're keeping
    reduced_f2_source_components = [
        f2_source_components[i] for i in range(len(f2_source_components)) if f2_name_components[i] not in f1_name_components
    ]
    source = '::'.join(f1_source_components + reduced_f2_source_components)

    if 'product' in feature1.qualifiers:
        f1_product_components = feature1.qualifiers['product'][0].split(' / ')
    else:
        f1_product_components = []
    if 'product' in feature2.qualifiers:
        f2_product_components = feature2.qualifiers['product'][0].split(' / ')
    else:
        f2_product_components = []
    reduced_f2_product_components = [_ for _ in f2_product_components if _ not in f1_product_components]
    product = ' / '.join(f1_product_components + reduced_f2_product_components)

    return name, source, product

def fusion_upgrade(base, upstream, downstream, update_location=False):
    """
    Update a feature to represent a fusion. Modifies the first argument.
    :param base: AutarkicSeqFeature the feature to change to represent the fusion.
    :param upstream: AutarkicSeqFeature the upstream component (left-most if + strand, right-most if - strand)
    :param downstream: gene: AutarkicSeqFeature
    :param update_location: bool whether to update the location to the outer termini of upstream and downstream.
    :returns: AutarkicSeqFeature base
    """
    gene, source, product = fusion_name(upstream, downstream)
    base.qualifiers['gene'][0] = gene
    base.source = source
    if product:
        base.qualifiers['product'][0] = product

    if update_location:
        if base.location.strand == 1:
            new_start = upstream.location.parts[0].start
            new_end = downstream.location.parts[-1].end
        else:
            new_start = downstream.location.parts[0].start
            new_end = upstream.location.parts[-1].end
        update_termini(base.location, new_start, new_end)

    base.bok = None
    pseudoscan.reset_pseudo(base)
    base.rcs, base.rce = upstream.rcs, downstream.rce
    base.de = False # With the fusion identity, it's no longer a delayed stop
    pseudoscan.call(base, ref_was_pseudo=False)

    return base

def fusionfisher(feature_list, adjudicate=True):
    """
    This function parses through a list of CDSs and returns a unique list of CDSs, cleaning up annotation artifacts due to gene fusion events, as well as renaming such genes and those that have conjoined with their neighbor.

    :param feature_list: list of sorted SeqFeature objects.
    :param adjudicate: bool whether to attempt resolution of detected conflicts.
       This is useful to have as True during RATT postprocessing, but to set False and leave for the thunderdome otherwise.
    :return:
       - list of the SeqFeature objects that were identified and tagged.
       - list of fusion genes
       - list of tuples of rejected SeqFeatures and a string describing why
    """
    logger = logging.getLogger('FusionFisher')
    outlist = []
    last_feature_by_strand = {}
    remarkable = {
        'hybrid': [],
        'conjoined': [],
    }
    rejects = []

    for feature in feature_list:
        outlist.append(feature)
        if feature.type != 'CDS':
            continue
        if feature.location.strand in last_feature_by_strand:
            prev_feature = last_feature_by_strand[feature.location.strand]
            pf_goodstart, pf_goodstop = prev_feature.rcs, prev_feature.rce
            cf_goodstart, cf_goodstop = feature.rcs, feature.rce
            last_feature_by_strand[feature.location.strand] = feature
        else:
            last_feature_by_strand[feature.location.strand] = feature
            continue


        if prev_feature.location == feature.location:
            #
            # Artifact: same gene name or gene name already included in fusion name string
            #           (i.e.,  geneA::geneB vs. geneA. We don't need want to make it geneA::geneB::geneA.)
            #
            if set(extractor.get_gene(prev_feature).split('::')).intersection(
                    set(extractor.get_gene(feature).split('::'))):
                if len(extractor.get_gene(feature)) > len(extractor.get_gene(prev_feature)):
                    goner = prev_feature
                    keeper = feature
                    outlist.remove(prev_feature)
                else:
                    goner = outlist.pop()
                    keeper = prev_feature
                rejects.append({
                    'feature':goner,
                    'superior':keeper,
                    'evid':'redundant_fusion_member',
                    'remark':"Same locus as rival (fusion) gene and name already included as fusion element there.",
                })
            #
            # Artifact due to a gene fusion hybrid
            #
            else:
                if (pf_goodstart and not cf_goodstart) or (not pf_goodstop and cf_goodstop):
                    upstream = prev_feature
                    downstream = feature
                elif (not pf_goodstart and cf_goodstart) or (pf_goodstop and not cf_goodstop):
                    upstream = feature
                    downstream = prev_feature
                else:
                    logger.warning(f"Could not determine which gene is first in fusion for {extractor.get_gene(prev_feature)} and {extractor.get_gene(feature)}. Using this order.")
                    upstream = prev_feature
                    downstream = feature

                fusion_upgrade(
                    base=prev_feature,
                    upstream=upstream,
                    downstream=downstream,
                )

                rejects.append({
                    'feature':outlist.pop(),
                    'superior':prev_feature,
                    'evid':'combined_annotation',
                    'remark':"Apparent hybrid fusion gene. Name incorporated into rival feature's and redundant locus removed.",
                })
                remarkable['hybrid'].append(prev_feature)

        elif have_same_stop(prev_feature.location, feature.location):
            #
            # Conjoined genes
            #
            if (((len(prev_feature.location) > len(feature.location)) and prev_feature.de)
                or ((len(feature.location) > len(prev_feature.location)) and feature.de)):
                if feature.location.strand == -1:
                    upstream = max(feature, prev_feature, key=lambda _: _.location.end)
                    downstream = min(feature, prev_feature, key=lambda _: _.location.end)
                else:
                    upstream = min(feature, prev_feature, key=lambda _: _.location.start)
                    downstream = max(feature, prev_feature, key=lambda _: _.location.start)

                designator.append_qualifier(
                    downstream.qualifiers,
                    'note',
                    f"Upstream gene {'|'.join([extractor.get_ltag(upstream),extractor.get_gene(upstream)])} conjoins with this one."
                )

                fusion_upgrade(
                    base=upstream,
                    upstream=upstream,
                    downstream=downstream,
                )

                remarkable['conjoined'].append(upstream)
            #
            # Another signature of a gene fusion hybrid
            #
            elif complementary(prev_feature, feature):
                if feature.rcs:
                    upstream, downstream = feature, prev_feature
                else:
                    upstream, downstream = prev_feature, feature
                fusion_upgrade(
                    base=upstream,
                    upstream=upstream,
                    downstream=downstream,
                    update_location=True,
                )

                remarkable['hybrid'].append(upstream)
                rejects.append({
                    'feature':downstream,
                    'superior':upstream,
                    'evid':'combined_annotation',
                    'remark':"Apparent hybrid fusion gene.",
                })
                outlist.remove(downstream)
            #
            # potential RATT misannotation, similar to the one in issue #51
            #
            elif adjudicate:
                if pf_goodstop and not cf_goodstop:
                    rejects.append({
                        'feature':outlist.pop(),
                        'superior':prev_feature,
                        'evid':'putative_misannotation',
                        'remark':"Has no reference-corresponding stop, while rival feature does, and both share the same stop position.",
                    })
                elif not pf_goodstop and cf_goodstop:
                    rejects.append({
                        'feature':prev_feature,
                        'superior':feature,
                        'evid':'putative_misannotation',
                        'remark':"Has no reference-corresponding stop, while rival feature does, and both share the same stop position.",
                    })
                    outlist.remove(prev_feature)
                #
                # likely scenarios in the case of a misannotation coinciding with a truncated gene
                #
                elif not pf_goodstart and cf_goodstart:
                    rejects.append({
                        'feature':prev_feature,
                        'superior':feature,
                        'evid':'putative_misannotation',
                        'remark':"Has no reference-corresponding coordinates, while rival feature has a reference-corresponding start, and both share the same stop position.",
                    })
                    outlist.remove(prev_feature)
                elif pf_goodstart and not cf_goodstart:
                    rejects.append({
                        'feature':outlist.pop(),
                        'superior':prev_feature,
                        'evid':'putative_misannotation',
                        'remark':"Has no reference-corresponding coordinates, while rival feature has a reference-corresponding start, and both share the same stop position."
                    })
                # remaining scenarios:
                # - both sets of coords are all good.
                # - both sets of coords are all bad.
                #      We could make a case for rejecting both, but that would seem to be creating a general rejection criterion which might not be warranted.
                #      I want to allow for the possibility of a double-truncation.
                # - both have good stops and bad starts
                # - both have bad stops and good starts
                #
                #
                # We will keep the longer feature unless only one of the two is pseudo, in which case we take the non-pseudo.
                else:
                    both_good_starts = pf_goodstart and cf_goodstart
                    both_good_stops = pf_goodstop and cf_goodstop
                    word_choice = lambda _: 'have' if _ else 'lack'
                    coord_status_report = (
                        f"Both {word_choice(both_good_starts)} "
                        f"reference-corresponding start codons and {word_choice(both_good_stops)} "
                        f"reference-corresponding stop codons.")
                    if designator.is_pseudo(prev_feature.qualifiers) != designator.is_pseudo(feature.qualifiers):
                        if designator.is_pseudo(prev_feature.qualifiers):
                            goner = prev_feature
                            keeper = feature
                            outlist.remove(prev_feature)
                        else:
                            goner = outlist.pop()
                            keeper = prev_feature
                        reason = "Non-pseudo"
                    elif len(prev_feature) > len(feature):
                        goner = outlist.pop()
                        keeper = prev_feature
                    else:
                        goner = prev_feature
                        keeper = feature
                        outlist.remove(prev_feature)
                    rejects.append({
                        'feature':goner,
                        'superior':keeper,
                        'evid':'shorter',
                        'remark':coord_status_report,
                    })

        if outlist[-1] != feature:
            last_feature_by_strand[feature.location.strand] = prev_feature


    return outlist, remarkable['hybrid'] + remarkable['conjoined'], rejects


def liftover_annotation(feature, ref_feature, inference):
    """
    Add ref_feature's functional annotation to feature.

    :param feature: SeqFeature ab initio annotation.
                    This argument is modified by this function.
    :param ref_feature: SeqFeature reference annotation
    :param inference: str /inference annotation justifying the liftover
    """

    for rubbish in ['gene', 'protein_id']:
        feature.qualifiers.pop(rubbish, None)
    # Remove inferences for the assignments that we're going to discard.
    # Only keep the ab initio inference from the ORF finder and ones that we assign.
    feature.qualifiers['inference'][:] = [
        _ for _ in feature.qualifiers['inference']
        if ('ab initio prediction' in _ or
            'Hybran' in _
            )
    ]
    # Add our own qualifier
    feature.qualifiers['inference'].append(
        inference
    )

    # make a copy of the reference qualifiers dict so
    # we can modify it before merging (we don't want to carry
    # over certain attributes)
    ref_feature_qualifiers_copy = deepcopy(ref_feature.qualifiers)
    ref_specific = [
        'locus_tag',
        'old_locus_tag',
        'translation',
        'pseudo',
        'pseudogene',
    ]
    for qual in ref_specific:
        ref_feature_qualifiers_copy.pop(qual, None)

    feature.qualifiers = merge_qualifiers(
        feature.qualifiers,
        ref_feature_qualifiers_copy,
    )

def thunderdome(abinit_annotation, ratt_annotation):
    """
    Two genes enter... one gene leaves.
    This function performs a 'standardized' comparison between RATT and Prokka
    annotations based on reference-correspondence and pseudo status. Will be used
    in check_inclusion_criteria for the cases dealing with conflicting annotations.
    :param abinit_annotation:
    :param ratt_annotation:
    :returns:
        - include_abinit (:py:class:`bool`) - whether the ab initio annotation should be kept
        - include_ratt (:py:class:`bool`) - whether the RATT annotation should be kept
        - evid (:py:class:`str`) - evidence code for rejected annotation
        - remark (:py:class:`str`) - explanation why the rejected annotation was not included
    """
    logger = logging.getLogger('Thunderdome')
    abinit_delayed_stop = abinit_annotation.de
    ratt_delayed_stop = ratt_annotation.de

    abinit_broken_stop = has_broken_stop(abinit_annotation)[0]
    ratt_broken_stop = has_broken_stop(ratt_annotation)[0]

    abinit_coord_status = coord_check(abinit_annotation, ref_annotation[
        key_ref_gene(abinit_annotation.source, abinit_annotation.qualifiers['gene'][0])
    ])
    ratt_coord_status = coord_check(ratt_annotation, ref_annotation[
        key_ref_gene(ratt_annotation.source, ratt_annotation.qualifiers['gene'][0])
    ])

    (abinit_start_ok, abinit_stop_ok) = abinit_coord_status
    abinit_coord_score = sum([int(_) for _ in abinit_coord_status])
    (ratt_start_ok, ratt_stop_ok) = ratt_coord_status
    ratt_coord_score = sum([int(_) for _ in ratt_coord_status])

    abinit_is_pseudo = designator.is_pseudo(abinit_annotation.qualifiers)
    ratt_is_pseudo = designator.is_pseudo(ratt_annotation.qualifiers)

    ratt_longer = (len(abinit_annotation.location) < len(ratt_annotation.location))
    abinit_longer = (len(abinit_annotation.location) > len(ratt_annotation.location))

    if abinit_is_pseudo == ratt_is_pseudo:
        # Both annotations being intact according to their respective reference names
        # suggests that the reference genes are highly similar.
        # RATT's assignment is furthermore based on synteny, so it wins out

        if ((all(ratt_coord_status) and all(abinit_coord_status)) or ratt_coord_status == abinit_coord_status):
            #Both non-pseudo and equal coord_status, use the more complete annotation.
            if not abinit_is_pseudo:
                include_abinit = False
                include_ratt = True
                evid = 'forfeit'
                remark = f"Equally valid call. RATT annotation is favored due to synteny."
            #Both pseudo
            else:
                if abinit_longer and abinit_delayed_stop and not ratt_delayed_stop:
                    include_abinit = True
                    include_ratt = False
                    evid = 'shorter_pseudo'
                    remark = "The ab initio annotation is favored over the RATT annotation because it is longer and contains a valid delayed stop."
                elif ratt_longer and ratt_delayed_stop and not abinit_delayed_stop:
                    include_abinit = False
                    include_ratt = True
                    evid = 'shorter_pseudo'
                    remark = "The RATT annotation is favored over the ab initio annotation because it is longer and contains a valid delayed stop."
                #Same location or equivalent delayed stop status
                else:
                    if abinit_broken_stop and not ratt_broken_stop:
                         include_abinit = False
                         include_ratt = True
                         evid = 'internal_stop'
                         remark = "The RATT annotation is favored over the ab initio annotation because it doesn't contain any internal stops."
                    elif ratt_broken_stop and not abinit_broken_stop:
                        include_abinit = True
                        include_ratt = False
                        evid = 'internal_stop'
                        remark = "The ab initio annotation is favored over the RATT annotation because it doesn't contain any internal stops."
                    else:
                        include_abinit = False
                        include_ratt = True
                        evid = 'forfeit'
                        remark = f"Equally valid call, but conflicts with RATT annotation; RATT favored due to synteny."

        elif (all(ratt_coord_status) or ratt_coord_score > abinit_coord_score or (ratt_stop_ok and not abinit_stop_ok)):
            include_abinit = False
            include_ratt = True
            evid = 'worse_ref_correspondence'
            remark = "RATT annotation more accurately named and delineated compared to the ab initio annotation."
        else:
        #This is the only other possibility:
        #elif (all(abinit_coord_status or abinit_coord_score > ratt_coord_score
        #or (abinit_stop_ok and not ratt_stop_ok)):
            include_abinit = True
            include_ratt = False
            evid = 'worse_ref_correspondence'
            remark = f"Ab initio annotation more accurately named and delineated compared to the RATT annotation."
    else:
        #Always take the non-pseudo annotation if possible
        if not abinit_is_pseudo and ratt_is_pseudo:
            include_abinit = True
            include_ratt = False
            evid = 'pseudo'
            remark = "Non-pseudo ab initio annotation takes precedence over the pseudo RATT annotation."
        else:
            include_abinit = False
            include_ratt = True
            evid = 'pseudo'
            remark = "Non-pseudo RATT annotation takes precedence over the pseudo ab initio annotation."

    if include_abinit == include_ratt:
        logger.warning(f"Both annotations were marked for {'inclusion' if include_ratt else 'exclusion'} and one annotation is expected to be excluded:\nRATT Feature\n{ratt_annotation}\n\nab initio Feature\n{abinit_annotation}\n\nUnhandled scenario--RATT favored due to synteny"
        )
        include_abinit = False
        include_ratt = True
        evid = 'forfeit'
        remark = "Both annotations marked for inclusion. Unhandled scenario--RATT annotation favored due to synteny"
    return include_abinit, include_ratt, evid, remark


def check_inclusion_criteria(
        ratt_annotation,
        abinit_annotation,
):
    """
    This function compares RATT and Prokka annotations and checks for conflicts.
    Either one feature or both will be accepted.
    If there is no conflict, both are kept. Otherwise, they are sent to the thunderdome().

    :param ratt_annotation:
    :param abinit_annotation:
    :returns:
        - include_abinit (:py:class:`bool`) - whether the ab initio annotation should be kept
        - include_ratt (:py:class:`bool`) - whether the RATT annotation should be kept
        - evid (:py:class:`str`) - evidence code for rejected annotation, if any. `None` otherwise.
        - remark (:py:class:`str`) - explanation why the rejected annotation, if any, was not included. `None` otherwise.
    """
    logger = logging.getLogger('CheckInclusionCriteria')
    include_ratt = True
    include_abinit = True
    evid = None
    remark = None

    if abinit_annotation.type != ratt_annotation.type:
        pass
    elif abinit_annotation.type != 'CDS':
        # TODO: we should come up with criteria for non-CDS genes
        if abinit_annotation.location == ratt_annotation.location:
            include_abinit = False
            evid = 'identical_non_cds'
    elif 'gene' not in abinit_annotation.qualifiers:
        if overlap_inframe(abinit_annotation.location, ratt_annotation.location):
            include_abinit = False
            evid = 'unnamed'
            remark = "Unnamed gene and conflicts (overlapping in-frame) with named rival annotation."
    else:
        same_gene_name = extractor.get_gene(ratt_annotation) == extractor.get_gene(abinit_annotation)
        same_loc = (abinit_annotation.location == ratt_annotation.location)
        if same_loc or same_gene_name:
            include_abinit, include_ratt, evid, remark = thunderdome(abinit_annotation, ratt_annotation)

        elif not same_gene_name and overlap_inframe(abinit_annotation.location, ratt_annotation.location):
            keepers, fusions, rejects = fusionfisher([ratt_annotation, abinit_annotation], adjudicate=False)
            for trial in rejects:
                reject = trial['feature']
                if reject == ratt_annotation:
                    include_ratt = False
                    evid = trial['evid']
                    remark = trial['remark']
                elif reject == abinit_annotation:
                    include_abinit = False
                    evid = trial['evid']
                    remark = trial['remark']
            # fusionfisher didn't detect a misannotation, but it didn't detect a fusion either.
            # welcome to the thunderdome!
            if not fusions and not rejects:
                include_abinit, include_ratt, evid, remark = thunderdome(abinit_annotation, ratt_annotation)

        else:
            #include everything if different names and not overlapping in frame
            include_abinit = True
            include_ratt = True

    return include_abinit, include_ratt, evid, remark

def fix_embl_id_line(embl_file):
    """

    :param embl_file:
    :return:
    """
    lines = []
    with open(embl_file, 'r') as embl:
        for line in embl:
            if line.startswith('ID'):
                lines.append(line.replace(' ; ; ; ; ; ', ' ; ; ; ; ; ; '))
            else:
                lines.append(line)
    with open(embl_file, 'w') as out:
        for line in lines:
            out.write(line)


def run(
        isolate_id,
        organism,
        strain,
        genome,
        annotation_fp,
        ref_proteins_fasta,
        ref_gbk_list,
        script_directory,
        seq_ident,
        seq_covg,
        ratt_enforce_thresholds,
        nproc=1,
):
    """
    Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file> -m
    <output_merged_genes> from the commandline. The log file output stats about the features that are added to the
    RATT annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/'isolate_id'.gbk and
    'isolate_id'/annomerge/'isolate_id'.log

    :param isolate_id: ID of the isolate (Example: H37Rv, 1-0006, etc.). This is the isolate_id that is used for naming
     Genbank files in Prokka
    :param organism: str binomial organism name or genus
    :param strain: str strain name
    :param genome: fasta file name corresponding to genome to annotate
    :param annotation_fp: Filepath where RATT, Prokka, reference and prokka no-reference annotations are located.
    Annomerge assumes that RATT annotations are located in <annotation_fp>/ratt, Prokka reference annotations are
    located in <annotation_fp>/prokka and prokka annotations without reference is located in
    <annotation_fp>/prokka-noreference. Additionally annomerge also assumes that withing prokka and prokka-noreference
    directories, the genbank files are located in <isolate_id>.gbk
    :param ref_proteins_fasta: File path for proteome fasta of reference strain
    :param ref_gbk_list: list of file paths for annotated GenBank file for reference genomes
    :param script_dir: Directory where hybran scripts are located
    :param ratt_enforce_thresholds: boolean - whether to enforce seq_ident/seq_covg for RATT-transferred annotations
    :param nproc: int number of processers available for use
    :return: EMBL record (SeqRecord) of annotated isolate
    """

    start_time = time.time()

    # avoid circular imports
    from . import ratt, prokka

    hybran_tmp_dir = config.hybran_tmp_dir
    global script_dir
    script_dir = script_directory
    global genetic_code
    genetic_code = config.cnf.genetic_code
    logger = logging.getLogger('Annomerge')

    if annotation_fp.endswith('/'):
        file_path = annotation_fp + isolate_id + '/'
    else:
        file_path = annotation_fp + '/' + isolate_id + '/'

    output_merged_genes = os.path.join(isolate_id, 'annomerge', 'merged_genes.gbk')
    output_genbank = os.path.join(isolate_id, 'annomerge', isolate_id + '.gbk')
    ratt_rejects = []
    ratt_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'ratt_unused.tsv')
    prokka_rejects = []
    prokka_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'prokka_unused.tsv')

    # create a dictionary of reference CDS annotations (needed for liftover to ab initio)
    global ref_annotation
    ref_annotation = keydefaultdict(ref_fuse)
    for ref_gbk_fp in ref_gbk_list:
        ref_id = os.path.basename(os.path.splitext(ref_gbk_fp)[0])
        for ref_record in SeqIO.parse(ref_gbk_fp, 'genbank'):
            ref_contig_id = '.'.join([ref_id, ref_record.id])
            for feature in ref_record.features:
                get_and_remove_ref_tracer(feature) # prevent our tracer note from propagating to future liftovers
                if feature.type != "CDS":
                    continue
                # setting feature.ref doesn't work for CompoundLocations
                for part in feature.location.parts:
                    part.ref = ref_contig_id
                feature.references = {ref_contig_id: ref_record.seq}
                # if reference paralogs have been collapsed, the last occurrence in the genome
                # will prevail.
                ref_annotation[key_ref_gene(ref_contig_id, feature.qualifiers['gene'][0])] = feature

    annomerge_records = list(SeqIO.parse(genome, "fasta"))
    contigs = [record.id for record in annomerge_records]

    ratt_file_path = os.path.join(file_path, 'ratt')
    ratt_features = ratt.postprocess(
        isolate_id,
        contigs,
        ratt_outdir=ratt_file_path,
        postprocess_outdir=os.path.join(file_path, 'ratt-postprocessed'),
        ref_annotation=ref_annotation,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
        nproc=nproc,
        enforce_thresholds=ratt_enforce_thresholds,
    )

    abinit_file_path = os.path.join(file_path, 'prokka')
    abinit_features = prokka.postprocess(
        isolate_id,
        contigs,
        prokka_outdir=abinit_file_path,
        postprocess_outdir=os.path.join(file_path, 'prokka-postprocessed'),
        ref_annotation=ref_annotation,
        ref_proteome=ref_proteins_fasta,
        seq_ident=seq_ident,
        seq_covg=seq_covg,
        nproc=nproc,
    )

    logger.info('Running Annomerge on ' + isolate_id)

    output_isolate_recs = []

    for i, contig in enumerate(contigs):
        seqname = '.'.join([isolate_id, contig])
        annomerge_records[i].id = annomerge_records[i].name = seqname
        annomerge_records[i].annotations['source'] = f"{organism} {strain}" if strain else organism
        annomerge_records[i].annotations['organism'] = organism
        annomerge_records[i].annotations['comment'] = "Annotated using hybran " + __version__ + " from https://lpcdrp.gitlab.io/hybran."
        annomerge_records[i].annotations['molecule_type'] = 'DNA'

        ratt_contig_features = ratt_features[contig]
        abinit_contig_features = abinit_features[contig]

        if len(ratt_contig_features) == 0:
            logger.info(f"{seqname}: Using ab initio annotations only since RATT did not annotate any")
            annomerge_records[i].features = abinit_contig_features
        elif len(abinit_contig_features) == 0:
            logger.info(f"{seqname}: Using RATT annotations only since ab initio methods did not annotate any")
            annomerge_records[i].features = ratt_contig_features
        else:
            annomerge_contig_features = []

            ratt_contig_features_dict = generate_feature_dictionary(ratt_contig_features)
            abinit_contig_features_dict = generate_feature_dictionary(abinit_contig_features)

            (
                colocated,
                inframe_conflicts,
                potential_unique_ratt,
                potential_unique_abinit,
                unique_ratt_features,
                unique_abinit_features,
                overlap_G,
            ) = compare(ratt_contig_features, abinit_contig_features)

            # Incorporate unique annotations categorically.
            annomerge_contig_features += unique_abinit_features
            logger.info(f"{seqname}: {len(unique_abinit_features)} ab initio features fall squarely into RATT's unannotated regions.")

            abinit_conflicts = collections.defaultdict(list)

            #
            # Redundant annotations
            #
            abinit_duplicates = []
            for ratt_feature, abinit_feature in colocated:
                if extractor.get_gene(ratt_feature) == extractor.get_gene(abinit_feature):
                    abinit_duplicates.append({
                        'feature': abinit_feature,
                        'superior': ratt_feature,
                        'evid': 'identical',
                    })
                else:
                    # We could call check_inclusion_criteria here immediately, but let's not change too much at once
                    abinit_conflicts[loc_index(abinit_feature)].append(loc_index(ratt_feature))
            prokka_rejects += abinit_duplicates
            logger.info(f"{seqname}: {len(abinit_duplicates)} ab initio features identical to RATT's")

            #
            # Overlapping, non-colocated annotations (both inframe and out-of-frame)
            #
            abinit_overlapping = potential_unique_abinit + [
                abinit_feature for (ratt_feature, abinit_feature) in inframe_conflicts
            ]
            for abinit_feature in abinit_overlapping:
                for _, ratt_feature_label in overlap_G.edges(abinit_feature.label):
                    ratt_feature = overlap_G.nodes[ratt_feature_label]['annotation']
                    abinit_conflicts[loc_index(abinit_feature)].append(loc_index(ratt_feature))

            logger.info(f"{seqname}: {len(inframe_conflicts)} ab initio ORFs conflicting in-frame with RATT's")
            logger.info(f"{seqname}: {len(abinit_conflicts)} ab initio features in total overlap RATT features. Resolving...")


            for feature_position in abinit_conflicts.keys():
                abinit_feature = abinit_contig_features_dict[feature_position]
                # Conflict Resolution
                # TODO: We're using set() here because the ratt_conflict_locs sometimes appear multiple times.
                #       This causes check_inclusion_criteria to run multiple times on the same pair, which
                #       in the case of gene fusions, causes cascading of the fusion names.
                #       (i.e.,  rattA, abinitB => rattA::abinitB => rattA::abinitB::abinitB => ...)
                for ratt_conflict_loc in set(abinit_conflicts[feature_position]):
                    # if the RATT annotation got rejected at some point, its remaining conflicts are moot
                    if ratt_conflict_loc not in ratt_contig_features_dict.keys():
                        include_abinit = True
                        continue
                    ratt_feature = ratt_contig_features_dict[ratt_conflict_loc]
                    include_abinit, include_ratt, evid, remark = check_inclusion_criteria(
                        ratt_annotation=ratt_feature,
                        abinit_annotation=abinit_feature,
                    )
                    if not include_abinit:
                        prokka_rejects.append({
                            'feature':abinit_feature,
                            'superior':ratt_feature,
                            'evid':evid,
                            'remark':remark,
                        })
                        break
                    # TODO: explain why this is an elif rather than an independent else.
                    elif not include_ratt:
                        ratt_rejects.append({
                            'feature': ratt_contig_features_dict.pop(ratt_conflict_loc),
                            'superior': abinit_feature,
                            'evid': evid,
                            'remark':remark,
                        })
                # Add the abinit feature if it survived all the conflicts
                if include_abinit:
                    annomerge_contig_features.append(abinit_feature)

            for ratt_feature_append in ratt_contig_features_dict.values():
                annomerge_contig_features.append(ratt_feature_append)

            annomerge_records[i].features = annomerge_contig_features


        # Finalize annotation records for this contig
        raw_features_unflattened = annomerge_records[i].features[:]
        raw_features = []
        for f_type in raw_features_unflattened:
            if isinstance(f_type, Bio.SeqFeature.SeqFeature):
                raw_features.append(f_type)
            elif isinstance(f_type, list) and len(f_type) > 0:
                for sub_feature in f_type:
                    if isinstance(sub_feature, Bio.SeqFeature.SeqFeature):
                        raw_features.append(sub_feature)
            else:
                continue
        annomerge_records[i].features = raw_features
        logger.debug(f'{seqname}: final feature annotation verification')
        n_final_cdss = 0
        for feature in annomerge_records[i].features:
            if feature.type == 'CDS':
                n_final_cdss += 1
                # Adding translated sequences where missing
                if 'translation' not in feature.qualifiers and not designator.is_pseudo(feature.qualifiers):
                    feature_sequence = translate(
                        feature.extract(),
                        table=genetic_code,
                        cds=True,
                    )
                    feature.qualifiers['translation'] = [feature_sequence]
                elif designator.is_pseudo(feature.qualifiers):
                    feature.qualifiers.pop('translation', None)
                if 'gene' not in feature.qualifiers.keys():
                    feature.qualifiers['gene'] = feature.qualifiers['locus_tag']

        sorted_final = get_ordered_features(annomerge_records[i].features)
        annomerge_records[i].features = sorted_final

        logger.info(f'{seqname}: {n_final_cdss} CDSs annomerge')

    with open(ratt_rejects_logfile, 'w') as ratt_rejects_log:
        log_feature_fates(ratt_rejects, logfile=ratt_rejects_log)
    with open(prokka_rejects_logfile, 'w') as prokka_rejects_log:
        log_feature_fates(prokka_rejects, logfile=prokka_rejects_log)

    SeqIO.write(annomerge_records, output_genbank, 'genbank')

    annomerge_records_dict = {i: annomerge_records[i].features for i in range(len(annomerge_records))}
    corrected_orf_logfile = os.path.join(
        isolate_id,
        'annomerge',
        'coord_corrections.tsv'
    )
    with open(corrected_orf_logfile, 'w') as corr_log:
        log_coord_corrections(annomerge_records_dict, corr_log)

    pseudoscan_logfile = os.path.join(
        isolate_id,
        'annomerge',
        'pseudoscan_report.tsv'
    )
    with open(pseudoscan_logfile, 'w') as p_log:
        log_pseudos(annomerge_records_dict, p_log)

    logger.debug('postprocessing and annomerge run time: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes')
