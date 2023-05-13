__author__ = "Deepika Gunasekaran"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@sdsu.edu"
__status__ = "Development"

# Title: Merge annotation from Prokka for the positions which are not annotated by RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files in case of
# multiple contigs/chromosome annotations and a Genbank file (.gbk) file from Prokka annotation run with a reference and
# a Genbank file (.gbk) file from Prokka annotation run without a reference. The output is an EMBL file with annotation
# predominantly from RATT and the intergenic regions annotated by RATT are filled with Prokka. This script also
# generates a log file to indicate characteristics of the transferred features from Prokka.

import sys
import functools
import collections
import os
import tempfile
import pickle
import logging
import time
import re
from copy import deepcopy
from math import log, ceil

# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing
import Bio
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition, SeqFeature, BeforePosition, AfterPosition
from Bio import Align

from . import BLAST
from . import converter
from . import config
from . import designator
from . import extractor
from . import __version__


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
        ref_contig_id = marker_note.split(':')[1]

    return ref_contig_id

def overlap_inframe(loc1, loc2):
    """
    Say whether two FeatureLocations are overlapping and share the same reading frame.
    This is for the purpose of CDSs to determine whether they should be corresponding to the
    same genomic feature or whether one gene has lost its stop codon and merged into its neighbor.

    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features overlap in-frame
    """
    # overlap determination adapted from https://stackoverflow.com/a/2953979
    overlap = (min(loc1.end, loc2.end) - max(loc1.start, loc2.start))

    if loc1.strand == loc2.strand and overlap > 0:
        # pseudogenes may occupy multiple reading frames, so
        # check both the start-defined and stop-defined frames,
        # as well as the actual overlapping part.
        if (loc1.start == loc2.start
            or loc1.end == loc2.end
            or (
                ((loc1.start - loc2.start) % 3 == 0
                 or (loc1.end - loc2.end) % 3 == 0)
                and (overlap % 3 == 0))
        ):
            return True
    return False

def have_same_stop(loc1, loc2):
    """
    Say whether two FeatureLocations have the same stop position.
    This is a special case of overlap_inframe().
    :param loc1: FeatureLocation
    :param loc2: FeatureLocation
    :return: True if both features have the same stop position
    """
    if loc1.strand == loc2.strand:
        return (
            (loc1.strand == 1 and loc1.end == loc2.end)
            or
            (loc1.strand == -1 and loc1.start == loc2.start)
        )
    return False

def has_delayed_stop(feature):
    """
    A SeqFeature with a "delayed stop" is determined by pseudoscan.
    The criteria required for having a delayed stop include the following:
    1) Divisible by three
    2) No internal/missing stops (is_broken_stop() = False)
    3) Good start and bad stop (coord_check() = True, False)
    4) Length of feature is longer than the reference gene
    :param feature: A SeqFeature object
    :return: True if the feature has a "delayed stop"
    """
    delayed_stop = []
    if 'note' in feature.qualifiers:
        delayed_stop = [_ for _ in feature.qualifiers['note'] if "delayed stop" in _]
    return bool(delayed_stop)

def has_valid_start(feature):
    """
    Finds the first three base pairs in a SeqFeaure and determines if it
    has a valid start codon according to its corresponding codon table.
    :param feature: A SeqFeature object
    :return: True if valid start codon exists
    """
    start_codons = CodonTable.generic_by_id[genetic_code].start_codons
    feature_seq = str(feature.extract(record_sequence))[:3]

    return feature_seq in start_codons

def has_broken_stop(feature):
    """
    Finds the amount and location of internal stops.
    :param feature: A SeqFeature object
    """
    internal_stop = False
    note = ''
    translation = str(feature.extract(record_sequence).translate(to_stop=False, table=genetic_code))
    num_stop = [i for i,e in enumerate(translation) if e == "*"]
    num_internal_stop = [i for i,e in enumerate(translation) if e == "*" and i != (len(translation)-1)]
    if len(num_internal_stop) >= 1 or translation[-1] != "*":
        internal_stop = True
        if len(num_stop) == 0:
            note = f"No stop codons detected in the translated sequence"
        else:
            note = f"Internal stop detected in the following codon(s): {' '.join([str(i) for i in num_internal_stop])}"
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
    return_feature = deepcopy(feature)
    rec_len = len(record_sequence) -1

    if feature.strand == 1:
        extended_feature_start = feature.location.start
        extended_feature_end = len(record_sequence) - 1
    else:
        extended_feature_start = 0
        extended_feature_end = feature.location.end

    extended_feature.location = FeatureLocation(
        int(extended_feature_start),
        int(extended_feature_end),
        strand=feature.strand
    )
    if circularize and ((0 in feature) or (rec_len in feature)):
        circ = FeatureLocation(int(0), int(rec_len), strand=feature.strand)
        #Creates CompoundLocation
        if feature.strand == 1:
            extended_feature.location = extended_feature.location + circularize
        else:
            extended_feature.location = circularize + extended_feature.location
    extended_seq = extended_feature.extract(record_sequence).translate(to_stop=True, table=genetic_code)

    #Translation extends up to, but not including the stop codon, so we need to add 3 to the end
    if feature.strand == 1:
        extended_feature_start = feature.location.start
        extended_feature_end = feature.location.start + (3*len(extended_seq)) + 3
        if extended_feature_end > rec_len:
            #Annotation will look like [feature.start ... >(rec_len)]
            extended_feature_end = AfterPosition(rec_len)
    else:
        extended_feature_start = feature.location.end - (3*len(extended_seq)) - 3
        extended_feature_end = feature.location.end
        if extended_feature_start < 0:
            #Annotation will look like [<0 ... feature.end]
            extended_feature_start = BeforePosition(0)

    return_feature.location = FeatureLocation(
        int(extended_feature_start),
        int(extended_feature_end),
        strand=feature.strand
    )
    return return_feature

def log_feature_fate(feature, logfile, remark=""):
    """
    General-purpose logging function to print out a gene's information and a comment
    :param feature: A SeqFeature object
    :param logfile: An open filehandle
    :param remark: (str) A comment
    """
    if 'locus_tag' in feature.qualifiers:
        locus_tag = feature.qualifiers['locus_tag'][0]
    else:
        locus_tag = feature.id
    print('\t'.join([locus_tag, remark]), file=logfile)

def log_coord_correction(og_feature, feature, logfile):
    """
    This function is used to log gene information when coord_check() fixes a start/stop postition
    :param og_feature: Original SeqFeature object
    :param feature: Updated SeqFeature object
    :param logfile: An open filehandle
    """
    locus_tag = og_feature.qualifiers['locus_tag'][0]
    gene_name = og_feature.qualifiers['gene'][0]
    strand = str(og_feature.strand)
    og_start = str(int(og_feature.location.start) + 1)
    og_end = (og_feature.location.end)
    new_start = str(int(feature.location.start) + 1)
    new_end = (feature.location.end)
    start_fixed = str(og_start != new_start).lower()
    stop_fixed = str(og_end != new_end).lower()

    if designator.is_pseudo(og_feature.qualifiers) == designator.is_pseudo(feature.qualifiers):
        if designator.is_pseudo(og_feature.qualifiers):
            status = "remains_pseudo"
        else:
            status = "remains_non_pseudo"
    elif designator.is_pseudo(feature.qualifiers):
        status = "became_pseudo"
    else:
        status = "became_non_pseudo"

    if og_feature.strand == -1:
        start_fixed, stop_fixed = stop_fixed, start_fixed
    line = [locus_tag, gene_name, strand, og_start, og_end, new_start, new_end,  start_fixed, stop_fixed, status]
    print('\t'.join(str(v) for v in line), file=logfile)


def load_reference_info(proteome_fasta):
    """
    This function gets some reference data for annomerge.
    1. A list of genes from the reference protein fasta
    2. A list of locus tags from the reference protein fasta
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    reference_gene_list = []
    reference_locus_list = []

    ref_fasta_records = SeqIO.parse(proteome_fasta, 'fasta')
    for record in ref_fasta_records:
        gene_locus_name = record.id
        gene_locus = gene_locus_name.split(':')
        if len(gene_locus[0]) == 0:  # If locus tag is not specified
            gene = gene_locus[1]
            locus = gene
        elif len(gene_locus[1]) == 0:  # If gene name is not specified
            locus = gene_locus[0]
            gene = locus
        else:
            locus = gene_locus[0]
            gene = gene_locus[1]
        reference_gene_list.append(gene)
        reference_locus_list.append(locus)
    return reference_gene_list, reference_locus_list


# Thanks to Jochen Ritzel
# https://stackoverflow.com/a/2912455
class keydefaultdict(collections.defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

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
        raise KeyError(fusion_gene_name)
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

def fissionfuser(flist, seq_ident, seq_covg, abinit_blast_results):
    """
    Given a list of features ordered by genomic position, identify adjacent gene fragments and combine them into a single feature.
    :param flist: list of SeqFeature objects
    :param seq_ident: sequence identity threshold for BLAST (for pseudo-calling)
    :param seq_covg: alignment coverage threshold for BLAST (for pseudo-calling)
    :param abinit_blast_results: dictionary of ab initio annotation locus tag to top blast hit stats (for pseudo-calling)
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
            if ((any(lg_status) and any(cg_status))
                and (int(lg_status[0])+int(cg_status[0]), int(lg_status[1])+int(cg_status[1]))==(1,1)):
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
            new_feature.location = FeatureLocation(new_start, new_end, feature.location.strand)
            new_feature.qualifiers = merge_qualifiers(dropped_feature.qualifiers, new_feature.qualifiers)

            if all(coord_check(
                    new_feature,
                    ref_annotation[key_ref_gene(new_feature.source, new_feature.qualifiers['gene'][0])],
                    fix_start=True, fix_stop=True)
                   ) or reason == 'overlapping_inframe':
                dropped_ltag_features.append(
                    (dropped_feature, f"{dropped_feature_name} combined with {new_feature_name}: {reason}")
                )
                # Re-call pseudoscan for updated notes and blast
                pseudoscan(
                    new_feature,
                    ref_annotation[key_ref_gene(new_feature.source, new_feature.qualifiers['gene'][0])],
                    seq_ident=seq_ident,
                    seq_covg=seq_covg,
                    blast_hit_dict=abinit_blast_results[new_feature.qualifiers['locus_tag'][0]],
                )
                last_gene_by_strand[feature.location.strand] = new_feature
                outlist.remove(last_gene)
                outlist.remove(feature)
                outlist.append(new_feature)
            else:
                logger.debug(f"Did not combine {dropped_feature_name} and {new_feature_name} ({reason}) due to failing coordinate verification.")

    return outlist, dropped_ltag_features


def fusionfisher(feature_list):
    """
    This function parses through a list of CDSs and returns a unique list of CDSs, cleaning up annotation artifacts due to gene fusion events, as well as renaming such genes and those that have conjoined with their neighbor.

    :param feature_list: list of sorted SeqFeature objects.
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
                rejects.append((
                    goner,
                    f'Redundant annotation with {extractor.get_ltag(keeper)}:{extractor.get_gene(keeper)}'
                ))
            #
            # Artifact due to a gene fusion hybrid
            #
            else:
                (pf_goodstart, pf_goodstop) = coord_check(
                    prev_feature,
                    ref_annotation[key_ref_gene(prev_feature.source, extractor.get_gene(prev_feature))],
                )
                (cf_goodstart, cf_goodstop) = coord_check(
                    feature,
                    ref_annotation[key_ref_gene(feature.source, extractor.get_gene(feature))],
                )
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

                prev_feature.qualifiers['gene'][0] = '::'.join([extractor.get_gene(_) for _ in [upstream, downstream]])
                prev_feature.source = '::'.join([upstream.source, downstream.source])
                prev_feature.qualifiers['product'] = [' / '.join([_.qualifiers['product'][0] for _ in [upstream, downstream] if 'product' in _.qualifiers])]

                rejects.append((
                    outlist.pop(),
                    f"Apparent hybrid fusion gene. Removed due to redundant location with {extractor.get_ltag(prev_feature)}:{extractor.get_gene(prev_feature)}.",
                ))
                remarkable['hybrid'].append(prev_feature)

        elif have_same_stop(prev_feature.location, feature.location):
            #
            # Conjoined genes
            #
            if (((len(prev_feature.location) > len(feature.location)) and has_delayed_stop(prev_feature))
                or ((len(feature.location) > len(prev_feature.location)) and has_delayed_stop(feature))):
                if feature.strand == -1:
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

                upstream.qualifiers['gene'][0] = '::'.join([extractor.get_gene(_) for _ in [upstream, downstream]])
                upstream.source = '::'.join([upstream.source, downstream.source])
                upstream.qualifiers['product'] = [' / '.join([_.qualifiers['product'][0] for _ in [upstream, downstream] if 'product' in _.qualifiers])]

                remarkable['conjoined'].append(upstream)
            #
            # potential RATT misannotation, similar to the one in issue #51
            #
            else:
                (pf_goodstart, pf_goodstop) = coord_check(
                    prev_feature,
                    ref_annotation[
                        key_ref_gene(prev_feature.source, extractor.get_gene(prev_feature))
                    ],
                )
                (cf_goodstart, cf_goodstop) = coord_check(
                    feature,
                    ref_annotation[
                        key_ref_gene(feature.source, extractor.get_gene(feature))
                    ],
                )
                if pf_goodstop and not cf_goodstop:
                    rejects.append((
                        outlist.pop(),
                        f"putative misannotation: has no reference-corresponding stop, while {extractor.get_ltag(prev_feature)}:{extractor.get_gene(prev_feature)} does, and both share the same stop position."
                    ))
                elif not pf_goodstop and cf_goodstop:
                    rejects.append((
                        prev_feature,
                        f"putative misannotation: has no reference-corresponding stop, while {extractor.get_ltag(feature)}:{extractor.get_gene(feature)} does, and both share the same stop position."
                    ))
                    outlist.remove(prev_feature)
                #
                # likely scenarios in the case of a misannotation coinciding with a truncated gene
                #
                elif not pf_goodstart and cf_goodstart:
                    rejects.append((
                        prev_feature,
                        f"putative misannotation: has no reference-corresponding coordinates, while {extractor.get_ltag(feature)}:{extractor.get_gene(feature)} has a reference-corresponding start, and both share the same stop position."
                    ))
                    outlist.remove(prev_feature)
                elif pf_goodstart and not cf_goodstart:
                    rejects.append((
                        outlist.pop(),
                        f"putative misannotation: has no reference-corresponding coordinates, while {extractor.get_ltag(prev_feature)}:{extractor.get_gene(prev_feature)} has a reference-corresponding start, and  both share the same stop position."
                    ))
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
                        f"Both {extractor.get_ltag(prev_feature)}:{extractor.get_gene(prev_feature)} and "
                        f"{extractor.get_ltag(feature)}:{extractor.get_gene(feature)} {word_choice(both_good_starts)} "
                        f"reference-corresponding start codons and {word_choice(both_good_stops)} "
                        f"reference-corresponding stop codons.")
                    reason = 'Longer'
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
                    rejects.append((
                        goner,
                        f'{coord_status_report} {reason} feature {extractor.get_ltag(keeper)}:{extractor.get_gene(keeper)} favored.'
                    ))

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

def coord_check(feature, ref_feature, fix_start=False, fix_stop=False, ref_gene_name=None
):
    """
    This function takes a feature as an input and aligns it to the corresponding reference gene.
    This function will return two Boolean values (True if the start/end coordinates align with the reference).
    The start/end coordinates can be corrected to match the reference if desired in various circumstances.
    :param feature: SeqFeature object
    :param fix_start: Boolean
    :param fix_stop: Boolean
    :param ref_gene_name: str reference gene to check against.
        Must be a key existing in the `ref_annotation` dictionary.
        If not defined, the reference gene matching `feature`'s gene qualifier is used instead.
    :return: True/False if the start/stop was fixed
    """

    logger = logging.getLogger('CoordCheck')
    ref_seq = extractor.get_seq(ref_feature)
    ref_length = len(ref_seq)
    feature_start = int(feature.location.start)
    feature_end = int(feature.location.end)
    feature_seq = feature.extract(record_sequence)
    og_feature = deepcopy(feature)
    if 'gene' not in og_feature.qualifiers:
        og_feature.qualifiers['gene'] = [ref_gene_name]
    og_feature_start = int(og_feature.location.start)
    og_feature_end = int(og_feature.location.end)

    def coord_align(ref_seq, feature_seq):
        #Probability that the continuous interval used to find good start/stops
        #occurs by chance should be = 1/(ref_seq*10)
        interval = max(ceil((log(len(ref_seq) * 10))/log(4)), 3)

        #The pairwise aligner returns an iterable list of alignments all with the
        #highest score it could possibly produce. Take the first 5 alignments from this list.
        aligner = Align.PairwiseAligner(scoring="blastn", mode = 'global')
        alignment = aligner.align(ref_seq, feature_seq)
        alignment_list = []
        while len(alignment_list) < 5:
            i = len(alignment_list)
            try:
                target = alignment[i].aligned[0]
                query = alignment[i].aligned[1]
                found_low = (target[0][0] == 0) and (abs(target[0][0] - target[0][1])) >= interval
                found_high = (target[-1][1] == ref_length) and (abs(target[-1][0] - target[-1][1])) >= interval
                alignment_list.append([i, len(target), len(query), found_low, found_high])
                continue
            except IndexError:
                break

        #If possible, filter to find the best fitting alignment rather than the first one available.
        if len(alignment_list) > 1:
            max_found = max(map(lambda x: [x[3],x[4]], alignment_list))
            alignment_list = [i for i in alignment_list if sum([i[3], i[4]]) == sum(max_found)]
            min_breaks = min(map(lambda x: [x[1],x[2]], alignment_list))
            alignment_list = [i for i in alignment_list if sum([i[1], i[2]]) == sum(min_breaks)]
            alignment = alignment[alignment_list[0][0]]
        else:
            alignment = alignment[0]

        target = alignment.aligned[0]
        query = alignment.aligned[1]
        score = alignment.score
        padding = False

        found_low = (target[0][0] == 0) and (abs(target[0][0] - target[0][1])) >= interval
        found_high = (target[-1][1] == ref_length) and (abs(target[-1][0] - target[-1][1])) >= interval

        target_low_seq = ref_seq[target[0][0]:target[0][1]]
        target_high_seq = ref_seq[target[-1][0]:target[-1][1]]
        query_low_seq = feature_seq[query[0][0]:query[0][1]]
        query_high_seq = feature_seq[query[-1][0]:query[-1][1]]

        #make sure there aren't too many mismatches causing falsely assigned found_low/high values
        if (target_low_seq[:3] != query_low_seq[:3]):
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
        return found_low, found_high, target, query, alignment, padding, score, interval

    def add_padding(feature, target, query, interval):
        #If we're looking to make corrections, add some context to help the aligner
        pad_feature = deepcopy(feature)
        feature_start = feature.location.start
        feature_end = feature.location.end

        #The only time padding matters is when the reference overhangs on the left or right side of
        #the feature. All other scenarios would be where found_low/high = True and extra context is unnecessary.
        left_ref_overhang = int(target[0][0]) > int(query[0][0])
        right_ref_overhang = (int(query[-1][1]) == len(feature)) and (int(target[-1][1]) < len(ref_seq))
        left_pad = right_pad = abs(len(ref_seq) - len(og_feature.extract(record_sequence))) + 2*interval

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

        pad_feature.location = FeatureLocation(
            feature_start,
            feature_end,
            strand=feature.strand
        )
        return pad_feature

    #First alignment
    found_low, found_high, target, query, alignment, padding, first_score, interval = coord_align(ref_seq, feature_seq)
    corrected_feature = deepcopy(feature)
    corrected_feature_start = corrected_feature.location.start
    corrected_feature_end = corrected_feature.location.end
    if padding:
        #Align again after adding padding to the feature sequence if warranted
        pad_feature = add_padding(feature, target, query, interval)
        pad_feature_seq = pad_feature.extract(record_sequence)

        pad_found_low, pad_found_high, pad_target, pad_query, pad_alignment, padding, second_score, second_interval = coord_align(ref_seq, pad_feature_seq)
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
            feature.location = FeatureLocation(
                int(og_feature_start),
                int(og_feature_end),
                strand=feature.strand
            )
            logger.warning(f"Attempted to correct {feature.qualifiers['gene'][0]} with invalid coordinates. Restoring original positions.")
        else:
            feature.location = FeatureLocation(
                int(feature_start),
                int(feature_end),
                strand=feature.strand
            )
        feature_seq = feature.extract(record_sequence)

        if i == 1:
            continue
        elif any([pad_found_low, pad_found_high]) and any([fix_start, fix_stop]) and (first_score > second_score):

            #Catch corner cases where we try to correct a reference corresponding start PAST the end of the feature.
            #or where we try to correct a reference corresponding stop BEFORE the start of a feature.
            #If the original alignment was really far off, and we found a downstream start/upstream stop
            #simply by chance from the addition of extra padding, we should just leave the gene alone.
            if corrected_feature_start > corrected_feature_end:
                 feature.location = FeatureLocation(
                     int(og_feature_start),
                     int(og_feature_end),
                     strand=feature.strand
                 )
                 logger.warning(f"Attempted to correct {feature.qualifiers['gene'][0]} with invalid coordinates. Restoring original positions.")
                 break
            corrected_feature.location = FeatureLocation(
                int(corrected_feature_start),
                int(corrected_feature_end),
                strand=feature.strand
            )
            corrected_feature_seq = corrected_feature.extract(record_sequence)
            cor_low, cor_high, cor_target, cor_query, cor_alignment, cor_padding, third_score, third_interval = coord_align(ref_seq, corrected_feature_seq)

            if (third_score >= first_score):
                second_score = third_score + 1
                continue
            else:
                break
        else:
            break

    #For cases where we DON'T want to fix coords, good_start/stop should only be true if positionally
    #identical to the start/stop in the reference. This is opposed to the normal definition of
    #good_start/stop as CONTAINING a reference corresponding start/stop site.
    if good_start and not fix_start:
        if ((feature.strand == 1 and feature.location.start != corrected_feature_start) or
            (feature.strand == -1 and feature.location.end != corrected_feature_end)):
            good_start = False
    if good_stop and not fix_stop:
        if ((feature.strand == 1 and feature.location.end != corrected_feature_end) or
            (feature.strand == -1 and feature.location.start != corrected_feature_start)):
            good_stop = False

    if og_feature.location != feature.location:
        feature.qualifiers['translation'] = [
            str(translate(
                feature.extract(record_sequence),
                table=genetic_code,
                to_stop=True,
            ))
        ]
        designator.append_qualifier(
            feature.qualifiers, 'inference',
            "COORDINATES:alignment:Hybran"
        )
        corrected_orf_report.append([og_feature, deepcopy(feature)])
    return good_start, good_stop

def pseudoscan(feature, ref_feature, seq_ident, seq_covg, attempt_rescue=False, blast_hit_dict=None
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
        previous_run_notes = [i for i in feature.qualifiers['note'] if 'Hybran/Pseudoscan' in i]
        # We're about to decide for ourselves whether the gene is pseudo in this current run
        # and we don't want the existence of the qualifier from the previous run
        # confounding the ref_was_pseudo determination.
        if previous_run_notes:
            feature.qualifiers['note'] = [i for i in feature.qualifiers['note'] if 'Hybran/Pseudoscan' not in i]
            feature.qualifiers.pop('pseudo', None)
            feature.qualifiers.pop('pseudogene', None)

    # checking the feature's pseudo attribute to decide whether the reference is pseudo
    # is only possible with RATT annotations since liftover from the reference may bring
    # the reference's pseudo tag attribute along with it. Our liftover doesn't do that.
    ref_was_pseudo = pseudo_note
    valid_start = has_valid_start(feature)
    og_broken_stop, og_stop_note = has_broken_stop(feature)
    divisible_by_three = lambda  _: len(_.location) % 3 == 0
    og_feature = deepcopy(feature)
    if ref_was_pseudo:
        good_start, good_stop = coord_check(feature, ref_feature)
        coords_ok = [good_start, good_stop]
        new_note = []
        new_note.append(f"Reference gene is pseudo")

        if all(coords_ok):
            new_note.append(f"Has reference-corresponding start and stop")

        if divisible_by_three(og_feature):
            if not divisible_by_three(ref_feature):
                new_note.append(f"Sequence divisible by three while the reference's is not.")
            else:
                new_note.append(f"Both this sequence and the reference's are divisible by three.")
        else:
            if not divisible_by_three(ref_feature):
                new_note.append(f"Both this sequence and the reference's have invalid reading frames -- not divisible by three.")
            else:
                new_note.append(f"Sequence not divisable by three while the reference sequence's is.")

        if (not all(coords_ok)) and (divisible_by_three(og_feature) and not divisible_by_three(ref_feature)) and not og_broken_stop:
            feature.qualifiers.pop('pseudo', None)
            feature.qualifiers.pop('pseudogene', None)
            is_pseudo = False
        else:
            if og_stop_note:
                new_note.append(f"{og_stop_note}")
            is_pseudo = True
            feature.qualifiers['pseudo'] = ['']

        new_note = ' | '.join(new_note)
        designator.append_qualifier(feature.qualifiers, 'note', f'Hybran/Pseudoscan: {new_note}')

    else:
        fix_start = False
        fix_stop = False
        confirmed_feature = False
        og_blast_defined = False
        while not confirmed_feature:
            good_start, good_stop = coord_check(feature, ref_feature, fix_start, fix_stop)
            coords_ok = [good_start, good_stop]
            broken_stop, stop_note = has_broken_stop(feature)

            if (feature.location != og_feature.location) and "No stop codons detected" in stop_note:
                extended_feature = stopseeker(feature)
                if len(extended_feature) > len(feature):
                    feature.location = extended_feature.location
                    good_start, good_stop = coord_check(feature, ref_feature)
                    coords_ok = [good_start, good_stop]
                    broken_stop, stop_note = has_broken_stop(feature)

            ref_seq = translate(
                extractor.get_seq(ref_feature),
                table=genetic_code, to_stop=True
            )
            feature_seq = translate(feature.extract(record_sequence), table=genetic_code, to_stop=True)
            #ref_match with 'thresholds enforced'
            top_hit, low_covg, blast_stats = BLAST.reference_match(
                query=SeqRecord(feature_seq),
                subject=SeqRecord(Seq(ref_seq), id=ref_feature.qualifiers['gene'][0]),
                seq_ident=seq_ident,
                seq_covg=seq_covg,
            )
            blast_ok = top_hit and not low_covg
            if og_blast_defined:
                lost_match = (og_blast_ok and not blast_ok)

            broken_stop, stop_note = has_broken_stop(feature)

            if not og_blast_defined:
                if all(coords_ok):
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
                # Revert coordinate correction if the blastp hit is lost.
                # (usually due to the reference-corresponding start being affected by an early frameshift)
                if lost_match:
                    feature.location = og_feature.location
                    feature.qualifiers = og_feature.qualifiers
                    broken_stop, stop_note = og_broken_stop, stop_note
                    blast_stats = og_blast_stats
                confirmed_feature = True

            if confirmed_feature:
                if blast_hit_dict:
                    blast_hit_dict.update(blast_stats[ref_feature.qualifiers['gene'][0]])
                # Summarize coord_check status for notes
                fancy_string = f"{'start' if not coords_ok[0] else ''}{' and ' if not any(coords_ok) else ''}{'end' if not coords_ok[1] else ''}"
                # Notes that are only interesting if we end up tagging the gene a certain way.
                new_note = []

                div_note = ('Locus has invalid reading frame-- not divisible by three', 0)
                if divisible_by_three(feature):
                    div_note = ('Locus divisible by three', 1)

                broke_note = ('No internal stop codons and ends with a valid stop codon', 0)
                if broken_stop:
                    broke_note = (stop_note, 1)

                start_note = ('Locus has invalid start codon', 0)
                if valid_start:
                    start_note = ('Locus has valid start codon', 1)

                coord_note = ('Locus has reference-corresponding start and end', list(map(int, coords_ok)))
                if not all(coords_ok):
                    coord_note = (f'Locus does not have reference-corresponding {fancy_string}', list(map(int, coords_ok)))

                blast_note = (f'Poor blastp match at {seq_ident}% identity and {seq_covg}% coverage thresholds', 0)
                if blast_ok:
                    blast_note = (f'Strong blastp match at {seq_ident}% identity and {seq_covg}% coverage thresholds', 1)

                #Note codes are categorized by a '0' or '1' and correspond to 'False' and 'True' respectively
                #D3 = Divisible by three [0/1]
                #VS = Valid start [0/1]
                #VE = Valid end [0/1]
                #RCS = Reference corresponding start [0/1]
                #RCE = Reference corresponding end [0/1]
                #BOK = Blast OK [0/1]
                note_codes = (f'D3{div_note[1]} VS{start_note[1]} VE{1 if broke_note[1] == 0 else 0} RCS{coord_note[1][0]} ' \
                              f'RCE{coord_note[1][1]} BOK{0 if not blast_ok else 1}')

                #This is the code for a normal non-pseudo gene, with no interesting characteristics.
                #These codes will not be added to the feature notes.
                if note_codes == 'D31 VS1 VE1 RCS1 RCE1 BOK1':
                    note_codes = None

                #The order in which notes appear is important. The first note should represent the
                #main reason for why a gene is marked pseudo or not.
                if ((all(coords_ok) or blast_ok) and divisible_by_three(feature)) and not broken_stop:
                    is_pseudo = False
                    if all(coords_ok):
                        if not blast_ok:
                            new_note.extend([coord_note[0], blast_note[0]])
                        if len(ref_feature) > len(feature):
                            new_note.append("Has deletion mutation(s) compared to the reference")
                            new_note.extend([broke_note[0], div_note[0]])
                        elif len(ref_feature) < len(feature):
                            new_note.append("Has insertion mutation(s) compared to the reference")
                            new_note.extend([broke_note[0], div_note[0]])
                    else:
                        #Strong blast, but non-reference-corresponding start/stop
                        new_note.extend([blast_note[0], coord_note[0]])
                else:
                    is_pseudo = True
                    feature.qualifiers['pseudo']=['']
                    #Primary reason for being pseudo. Followed by interesting cases if applicable.
                    if broken_stop or not divisible_by_three(feature): #Broken Stop/Not divisible by three
                        new_note.extend([broke_note[0], div_note[0]])

                        if all(coords_ok):
                            if not divisible_by_three(feature):
                                if len(ref_feature) > len(feature):
                                    new_note.append("Has deletion mutation(s) compared to the reference")
                                elif len(ref_feature) < len(feature):
                                    new_note.append("Has insertion mutation(s) compared to the reference")
                            new_note.append(coord_note[0])

                        if blast_ok:
                            new_note.append(blast_note[0])
                    #Valid reading frame AND not a broken stop.
                    #Must be pseudo because of (not all(coords_ok) AND not (blast_ok))
                    else:
                        if coords_ok == [True, False]:
                            if len(ref_feature) < len(feature):
                                new_note.append(coord_note[0])
                                new_note.append("Has a frameshift mutation leading to a delayed stop codon")
                        else:
                            new_note.append(coord_note[0])
                            new_note.append(blast_note[0])
                            new_note.append(broke_note[0])
                            new_note.append(div_note[0])

                if new_note:
                    new_note = ' | '.join(new_note)
                    designator.append_qualifier(feature.qualifiers, 'note', f'Hybran/Pseudoscan: {new_note}')
                if note_codes:
                    designator.append_qualifier(feature.qualifiers, 'note', note_codes)
                break
            else:
                og_blast_ok, og_blast_stats = (blast_ok, blast_stats)
                og_blast_defined = True

    return is_pseudo

def isolate_valid_ratt_annotations(feature_list, reference_locus_list, seq_ident, seq_covg, ratt_enforce_thresholds,
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

def find_inframe_overlaps(ratt_features, abinit_features_dictionary):
    """
    Identify the relevant ab initio features.
    Mark for removal those identically annotated by RATT and mark for conflict resolution those with different coordinates, but overlapping in frame.
    :param ratt_features: List of features from RATT (list of SeqFeature objects)
    :param abinit_features_dictionary: sorted dictionary of ab initio features ordered by the genomic position (i.e.,
    feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the corresponding
    SeqFeature object)
    :returns:
        - abinit_features_not_in_ratt (:py:class:`dict`) -
            sorted dictionary of ab initio features that are not annotated by RATT ordered by the genomic position
            i.e. feature location where key is a tuple (feature_start, feature_end, feature strand) and value is the
             corresponding SeqFeature object
        - ratt_overlapping_genes (:py:class:`dict`) -
            dictionary with ab initio feature location triples as keys and a list of RATT feature location triples
            as values. Conflicts here are annotations that overlap in-frame.
        - abinit_rejects (:py:class:`list`) -
            list of tuples of the form (abinit_feature, remark) where abinit_feature is a SeqFeature of a rejected
            ab initio annotation and remark is the rationale as a string.
    """

    abinit_features_not_in_ratt = abinit_features_dictionary.copy()
    ratt_overlapping_genes = collections.defaultdict(list)
    abinit_rejects = []
    for abinit_feature_position in abinit_features_dictionary.keys():
        abinit_feature = abinit_features_dictionary[abinit_feature_position]
        (abinit_start, abinit_end, abinit_strand) = abinit_feature_position
        abinit_duplicate_removed = False
        for ratt_feature in ratt_features:
            if ratt_feature.type != 'CDS':
                continue
            ratt_start = int(ratt_feature.location.start)
            ratt_end = int(ratt_feature.location.end)
            ratt_strand = ratt_feature.location.strand
            if ratt_start > abinit_end and ratt_end > abinit_end:
                break
            elif (overlap_inframe(ratt_feature.location, abinit_feature.location)):
                if (len(ratt_feature.location) == len(abinit_feature.location)
                    and extractor.get_gene(ratt_feature) == extractor.get_gene(abinit_feature)):
                    abinit_rejects.append((abinit_features_not_in_ratt.pop((abinit_start, abinit_end, abinit_strand), None),
                                           'duplicate of ' + ratt_feature.qualifiers['locus_tag'][0]))
                    abinit_duplicate_removed = True
                else:
                    ratt_overlapping_genes[abinit_feature_position].append((ratt_start, ratt_end, ratt_strand))
        if abinit_duplicate_removed:
            ratt_overlapping_genes.pop(abinit_feature_position, None)
    return abinit_features_not_in_ratt, ratt_overlapping_genes, abinit_rejects


def get_interregions(embl_record, intergene_length=1):
    """
    # Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
    # Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    # Do not remove this comment
    # This function was modified by Deepika Gunasekaran

    This function gets the genomic locations that do not have an coding-sequence (intergenic regions)
    :param embl_record: EMBL SeqRecord
    :param intergene_length: minimum length of integernic region (Default: 1)
    :return:
    SeqRecord of intergenic locations,
    list of intergenic positions where each element in the list is a tuple (start, end, strand),
    dictionary of genes preceding the intergenic regions where the key is a tuple (start of intergenic region, strand)
        and value is the SeqFeature preceding the intergenic region,
    dictionary of genes succeeding the intergenic regions where the key is a tuple (end of intergenic region, strand) and
        value is the SeqFeature succeeding the intergenic region
    """

    seq_record = embl_record
    cds_list_plus = []
    cds_list_minus = []
    intergenic_records = []
    intergenic_positions = []
    pre_intergene = {}
    post_intergene = {}
    # Loop over the genome file, get the CDS features on each of the strands
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        mystart = feature.location.start
        myend = feature.location.end
        if feature.strand == -1:
            cds_list_minus.append((mystart, myend, -1))
            pre_intergene[(myend, -1)] = feature
            post_intergene[(mystart, -1)] = feature
        elif feature.strand == 1:
            cds_list_plus.append((mystart, myend, 1))
            pre_intergene[(myend, 1)] = feature
            post_intergene[(mystart, 1)] = feature
        else:
            sys.stderr.write("No strand indicated %d-%d. Assuming +\n" % (mystart, myend))
            cds_list_plus.append((mystart, myend, 1))
            pre_intergene[(myend, 1)] = feature
            post_intergene[(mystart, 1)] = feature
    cds_list_plus = sorted(cds_list_plus)
    cds_list_minus = sorted(cds_list_minus)
    for i, pospair in enumerate(cds_list_plus[1:]):
        # Compare current start position to previous end position
        last_end = cds_list_plus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end: this_start]
            strand_string = "+"
            intergenic_records.append(SeqRecord(intergene_seq,
                                                id="%s-ign-%d" % (seq_record.name, i),
                                                description="%s %d-%d %s" % (seq_record.name, last_end + 1,
                                                                             this_start, strand_string)))
            intergenic_positions.append((last_end + 1, this_start, strand_string))
    for i, pospair in enumerate(cds_list_minus[1:]):
        last_end = cds_list_minus[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end: this_start]
            strand_string = "-"
            intergenic_records.append(SeqRecord(intergene_seq,
                                                id="%s-ign-%d" % (seq_record.name, i),
                                                description="%s %d-%d %s" % (seq_record.name, last_end + 1,
                                                                             this_start, strand_string)))
            intergenic_positions.append((last_end + 1, this_start, strand_string))
    return intergenic_records, intergenic_positions, pre_intergene, post_intergene


def populate_gaps(
        abinit_features,
        intergenic_positions,
        ratt_pre_intergene,
        ratt_post_intergene,
):
    """
    :param abinit_features: dict as from generate_feature_dictionary of ab initio features to try to incorporate
    :param intergenic_positions: (sorted) output of get_interregions
    """

    logger = logging.getLogger('PopulateGaps')
    logger.debug('Merging RATT and ab initio annotations')
    abinit_features = deepcopy(abinit_features)
    abinit_keepers = []
    abinit_conflicts = collections.defaultdict(list)
    prev_unannotated_region_end = {'+':1, '-':1}
    for j in intergenic_positions:
        # Variable definitions
        (ratt_unannotated_region_start, ratt_unannotated_region_end, ratt_strand) = j
        next_abinit_features = deepcopy(abinit_features)
        for feature_position in abinit_features.keys():
            prokka_feature = abinit_features[feature_position]
            (prokka_feature_start, prokka_feature_end, prokka_strand) = feature_position
            ratt_unannotated_region_range = range(ratt_unannotated_region_start,
                                                  ratt_unannotated_region_end + 1)
            if((prokka_strand == -1 and ratt_strand == '-')
               or (prokka_strand == 1 and ratt_strand == '+')):
                # If Prokka feature is location before the end of the previous intergenic region, pop the key from the
                # dictionary and continue loop
                if((prokka_feature_start < prev_unannotated_region_end[ratt_strand]
                    and prokka_feature_end <= prev_unannotated_region_end[ratt_strand])
                   or prokka_feature.type == 'source' # The source feature is accounted for by RATT
                ):
                    next_abinit_features.pop(feature_position, None)
                    continue
                # Else if the prokka feature location is after the end of the intergenic region, break out of
                # the inner loop
                elif(prokka_feature_start > ratt_unannotated_region_end
                     and prokka_feature_end > ratt_unannotated_region_end):
                    break
                # If the ab initio feature is contained in the unannotated range
                elif(prokka_feature_start in ratt_unannotated_region_range
                     and prokka_feature_end in ratt_unannotated_region_range):
                    abinit_keepers.append(prokka_feature)
                # If the Prokka feature overlaps with two RATT features
                elif(prokka_feature_start < ratt_unannotated_region_start
                     and prokka_feature_end > ratt_unannotated_region_end):
                    ratt_overlapping_feature_1 = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                     prokka_strand)]
                    ratt_overlapping_feature_2 = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                      prokka_strand)]
                    ratt_overlapping_feature_1_loc = (int(ratt_overlapping_feature_1.location.start),
                                                      int(ratt_overlapping_feature_1.location.end),
                                                      int(ratt_overlapping_feature_1.location.strand))
                    ratt_overlapping_feature_2_loc = (int(ratt_overlapping_feature_2.location.start),
                                                      int(ratt_overlapping_feature_2.location.end),
                                                      int(ratt_overlapping_feature_2.location.strand))
                    abinit_conflicts[feature_position] += [ratt_overlapping_feature_1_loc, ratt_overlapping_feature_2_loc]
                # If the Prokka feature overlaps with one RATT feature
                else:
                    does_not_overlap = False
                    # the abinit feature might overlap a ratt feature without crossing into the unannotated region.
                    # such cases are more likely to be conflicts, so we should identify them and get them resolved.
                    if (prokka_feature_start < ratt_unannotated_region_start) and \
                            (prokka_feature_end < ratt_unannotated_region_end):
                        ratt_overlapping_feature = ratt_pre_intergene[(ratt_unannotated_region_start - 1,
                                                                       prokka_strand)]
                    elif (prokka_feature_start in ratt_unannotated_region_range) and \
                            prokka_feature_end > ratt_unannotated_region_end:
                        ratt_overlapping_feature = ratt_post_intergene[(ratt_unannotated_region_end,
                                                                        prokka_strand)]
                    else:
                        does_not_overlap = True
                    if not does_not_overlap:
                        ratt_overlapping_feature_loc = (int(ratt_overlapping_feature.location.start),
                                                        int(ratt_overlapping_feature.location.end),
                                                        int(ratt_overlapping_feature.location.strand))
                        abinit_conflicts[feature_position].append(ratt_overlapping_feature_loc)
        prev_unannotated_region_end[ratt_strand] = ratt_unannotated_region_end
        abinit_features = next_abinit_features

    return abinit_keepers, abinit_conflicts

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
        - remark (:py:class:`str`) - explanation for why the rejected annotation, if any, was not included
    """
    logger = logging.getLogger('Thunderdome')
    abinit_delayed_stop = has_delayed_stop(abinit_annotation)
    ratt_delayed_stop = has_delayed_stop(ratt_annotation)

    abinit_broken_stop = has_broken_stop(abinit_annotation)[0]
    ratt_broken_stop = has_broken_stop(ratt_annotation)[0]

    abinit_coord_status = coord_check(abinit_annotation, ref_annotation[
        key_ref_gene(abinit_annotation.source, abinit_annotation.qualifiers['gene'][0])
    ])
    ratt_coord_status = coord_check(ratt_annotation, ref_annotation[
        key_ref_gene(ratt_annotation.source, ratt_annotation.qualifiers['gene'][0])
    ])

    #If an annotation has a delayed stop, (indicating a potential gene fusion event), we want to consider it a good stop
    #because it could correspond to the downstream gene's stop position. Gene fusions should take precedence over
    #alternative annotations if warranted.
    if abinit_delayed_stop:
        abinit_coord_status = (abinit_coord_status[0], True)
    if ratt_delayed_stop:
        ratt_coord_status = (ratt_coord_status[0], True)

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
                if abinit_longer:
                    include_abinit = True
                    include_ratt = False
                    remark = f"Equally valid call, but the more complete ab initio annotation is favored."
                elif ratt_longer:
                    include_abinit = False
                    include_ratt = True
                    remark = f"Equally valid call, but the more complete RATT annotation is favored."
                else:
                    include_abinit = False
                    include_ratt = True
                    remark = f"Equally valid call, but conflicts with RATT annotation {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)}; RATT favored due to synteny."
            #Both pseudo
            else:
                if abinit_longer and abinit_delayed_stop and not ratt_delayed_stop:
                    include_abinit = True
                    include_ratt = False
                    remark = f"The ab initio annotation is favored due to having a valid delayed stop."
                elif ratt_longer and ratt_delayed_stop and not abinit_delayed_stop:
                    include_abinit = False
                    include_ratt = True
                    remark = f"The ratt annotation is favored due to having a valid delayed stop."
                #Same location or equivalent delayed stop status
                else:
                    if abinit_broken_stop and not ratt_broken_stop:
                         include_abinit = False
                         include_ratt = True
                         remark = f"The ratt annotation is favored due to having a valid stop."
                    elif ratt_broken_stop and not abinit_broken_stop:
                        include_abinit = True
                        include_ratt = False
                        remark = f"The ab initio annotation is favored due to having a valid stop."
                    else:
                        include_abinit = False
                        include_ratt = True
                        remark = f"Equally valid call, but conflicts with RATT annotation {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)}; RATT favored due to synteny."

        elif (all(ratt_coord_status) or ratt_coord_score > abinit_coord_score or (ratt_stop_ok and not abinit_stop_ok)):
            include_abinit = False
            include_ratt = True
            remark = f"RATT annotation {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)} more accurately named and delineated."
        else:
        #This is the only other possibility:
        #elif (all(abinit_coord_status or abinit_coord_score > ratt_coord_score
        #or (abinit_stop_ok and not ratt_stop_ok)):
            include_abinit = True
            include_ratt = False
            remark = f"Ab initio annotation {extractor.get_ltag(abinit_annotation)}:{extractor.get_gene(abinit_annotation)} more accurately named and delineated."
    else:
        #Always take the non-pseudo annotation if possible
        if not abinit_is_pseudo and ratt_is_pseudo:
            include_abinit = True
            include_ratt = False
            remark = "Non-pseudo ab initio annotation takes precedence."
        else:
            include_abinit = False
            include_ratt = True
            remark = "Non-pseudo RATT annotation takes precedence."

    if include_abinit == include_ratt:
        logger.warning(f"Both annotations were marked for {'inclusion' if include_ratt else 'exclusion'} and one annotation is expected to be excluded:\nRATT Feature\n{ratt_annotation}\n\nab initio Feature\n{abinit_annotation}\n\nUnhandled scenario--RATT favored due to synteny"
        )
        include_abinit = False
        include_ratt = True
        remark = "Conflicting annotations, RATT favored due to synteny"
    return include_abinit, include_ratt, remark


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
        - remark (:py:class:`str`) - explanation for why the rejected annotation, if any, was not included
    """
    logger = logging.getLogger('CheckInclusionCriteria')
    include_ratt = True
    include_abinit = True
    remark = ''

    if abinit_annotation.type != ratt_annotation.type:
        pass
    elif abinit_annotation.type != 'CDS':
        # TODO: we should come up with criteria for non-CDS genes
        if abinit_annotation.location == ratt_annotation.location:
            include_abinit = False
            remark = f'Has same location as {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)}.'
    elif 'gene' not in abinit_annotation.qualifiers:
        if overlap_inframe(abinit_annotation.location, ratt_annotation.location):
            include_abinit = False
            remark = f"Hypothetical gene and conflicts (overlapping in-frame) with RATT's {extractor.get_ltag(ratt_annotation)}:{extractor.get_gene(ratt_annotation)}."
    else:
        same_gene_name = extractor.get_gene(ratt_annotation) == extractor.get_gene(abinit_annotation)
        same_loc = (abinit_annotation.location == ratt_annotation.location)
        if same_loc or same_gene_name:
            include_abinit, include_ratt, remark = thunderdome(abinit_annotation, ratt_annotation)

        elif not same_gene_name and overlap_inframe(abinit_annotation.location, ratt_annotation.location):
            keepers, fusions, rejects = fusionfisher([ratt_annotation, abinit_annotation])
            for (reject, reason) in rejects:
                if reject == ratt_annotation:
                    include_ratt = False
                    remark = reason
                elif reject == abinit_annotation:
                    include_abinit = False
                    remark = reason
            # fusionfisher didn't detect a misannotation, but it didn't detect a fusion either.
            # welcome to the thunderdome!
            if not fusions and not rejects:
                include_abinit, include_ratt, remark = thunderdome(abinit_annotation, ratt_annotation)

        else:
            #include everything if different names and not overlapping in frame
            include_abinit = True
            include_ratt = True

    return include_abinit, include_ratt, remark

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


def run(isolate_id, contigs, annotation_fp, ref_proteins_fasta, ref_gbk_list, script_directory, seq_ident, seq_covg, ratt_enforce_thresholds,
    nproc=1,
):
    """
    Annomerge takes as options -i <isolate_id> -g <output_genbank_file> -l <output_log_file> -m
    <output_merged_genes> from the commandline. The log file output stats about the features that are added to the
    RATT annotation. The default locations for the -g and -l options are 'isolate_id'/annomerge/'isolate_id'.gbk and
    'isolate_id'/annomerge/'isolate_id'.log

    :param isolate_id: ID of the isolate (Example: H37Rv, 1-0006, etc.). This is the isolate_id that is used for naming
     Genbank files in Prokka
    :param contigs: list of strings for the contig names
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

    hybran_tmp_dir = config.hybran_tmp_dir
    global script_dir
    script_dir = script_directory
    global genetic_code
    genetic_code = config.genetic_code
    logger = logging.getLogger('Annomerge')
    logger.debug('Running Annomerge on ' + isolate_id)
    start_time = time.time()

    reference_gene_list, reference_locus_list  = load_reference_info(ref_proteins_fasta)
    if annotation_fp.endswith('/'):
        file_path = annotation_fp + isolate_id + '/'
    else:
        file_path = annotation_fp + '/' + isolate_id + '/'
    ratt_file_path = file_path + 'ratt'
    ratt_correction_files = []
    ratt_gbk_files = {}
    try:
        for contig in contigs:
            embl_file = f"{isolate_id.replace('|','_')}.{contig.replace('|','_')}.final.embl"
            gbk = converter.convert_embl_to_gbk(ratt_file_path + '/' + embl_file)
            ratt_gbk_files[contig] = gbk
        correction_files = [cf for cf in os.listdir(ratt_file_path) if cf.endswith('.Report.txt')]
        for corr_file in correction_files:
            corr_file_path = ratt_file_path + '/' + corr_file
            ratt_correction_files.append(corr_file_path)
    except OSError:
        logger.error('Expecting RATT annotation files but found none')
    if not ratt_gbk_files:
        logger.error('RATT did not complete running. Please see the log for more details.')
    try:
        input_prokka_genbank = file_path + 'prokka/' + isolate_id + '.gbk'
    except OSError:
        logger.error('Expecting Prokka annotation file but found none')
    output_merged_genes = os.path.join(isolate_id, 'annomerge', 'merged_genes.gbk')
    output_genbank = os.path.join(isolate_id, 'annomerge', isolate_id + '.gbk')
    prokka_records = list(SeqIO.parse(input_prokka_genbank, 'genbank'))
    global ratt_rejects  # some RATT annotations are rejected as a side effect of check_inclusion_criteria
    ratt_rejects = []
    ratt_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'ratt_unused.tsv')
    prokka_rejects = []
    prokka_rejects_logfile = os.path.join(isolate_id, 'annomerge', 'prokka_unused.tsv')
    annomerge_records = []
    corrected_abinit_orf_logfile = os.path.join(isolate_id, 'prokka', 'hybran_coord_corrections.tsv')
    corrected_ratt_orf_logfile = os.path.join(isolate_id, 'ratt', 'hybran_coord_corrections.tsv')
    global corrected_orf_report
    corrected_orf_report = []
    # create a dictionary of reference CDS annotations (needed for liftover to ab initio)
    global ref_annotation
    ref_annotation = keydefaultdict(ref_fuse)
    for ref_gbk_fp in ref_gbk_list:
        ref_id = os.path.basename(os.path.splitext(ref_gbk_fp)[0])
        for ref_record in SeqIO.parse(ref_gbk_fp, 'genbank'):
            ref_contig_id = '.'.join([ref_id, ref_record.name])
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

    output_isolate_recs = []

    for i, contig in enumerate(contigs):
        seqname = '.'.join([isolate_id, contig])
        ratt_contig_record = SeqIO.read(ratt_gbk_files[contig], 'genbank')
        global record_sequence
        record_sequence = ratt_contig_record.seq
        prokka_contig_record = prokka_records[i]

        ratt_contig_features = ratt_contig_record.features
        prokka_contig_features = prokka_contig_record.features
        # When prokka assigns the same gene name to multiple orfs, it appends _1, _2, ... to make the names unique.
        # That causes issues for us because we expect all copies of a gene to have the same name.
        for f in prokka_contig_features:
            if 'gene' in f.qualifiers.keys():
                f.qualifiers['gene'][0] = re.sub(r"_\d+$","",f.qualifiers['gene'][0])


        #
        # RATT Postprocessing
        #
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
                                           reference_locus_list=reference_locus_list,
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

        if merged_features and i == 0:
            merged_features_record = prokka_contig_record[:]
            merged_features_record.features = merged_features
            SeqIO.write(merged_features_record, output_merged_genes, 'genbank')
        global ratt_contig_features_dict
        ratt_contig_features_dict = generate_feature_dictionary(ratt_contig_features)
        if len(ratt_contig_features) == 0:
            logger.warning(f"NO RATT ANNOTATION FOR {seqname}")
            feature_additions = {}
            feature_lengths = {}
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
            for prokka_feature in prokka_contig_record.features:
                if prokka_feature.type not in feature_additions.keys():
                    feature_additions[prokka_feature.type] = 1
                    feature_lengths[prokka_feature.type] = [len(prokka_feature.location)]
                else:
                    feature_additions[prokka_feature.type] += 1
                    feature_lengths[prokka_feature.type].append(len(prokka_feature.location))
            continue
        elif len(prokka_contig_features) == 0:
            logger.warning(f"NO AB INITIO ANNOTATION FOR {seqname}")
            prokka_contig_features = ratt_contig_features
            if len(merged_features) > 0:
                for feature in merged_features:
                    prokka_contig_features.append(feature)
            logger.warning(f'{seqname}: no ab initio annotations to add')
            prokka_contig_record.features = prokka_contig_features
            annomerge_records.append(prokka_contig_record)
        else:
            # Initializing annomerge gbf record to hold information such as id, etc from prokka but populating the
            # features from RATT
            add_prokka_contig_record = prokka_contig_record[:]
            add_prokka_contig_record.features = []

            try:
                ratt_contig_record_mod = ratt_contig_record[:]
            except AttributeError:
                logger.error('Contains features with fuzzy locations')
                logger.error(ratt_contig_record)
            ratt_contig_record_mod.features = ratt_contig_features

            #
            # Ab initio Postprocessing
            #
            abinit_features_dict = generate_feature_dictionary(prokka_contig_features)
            logger.debug(f'{seqname}: Checking ab initio CDS annotations for matches to reference using {nproc} process(es)')
            #
            # can contain results for hits to multiple reference genes
            abinit_blast_results_complete = {}
            # only contains results for the accepted reference gene hit
            abinit_blast_results = {}
            refmatch = functools.partial(
                BLAST.reference_match,
                subject=ref_proteins_fasta,
                seq_ident=seq_ident,
                seq_covg=seq_covg,
            )
            prokka_contig_cdss = [f for f in abinit_features_dict.values() if f.type == 'CDS']
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
                    ref_id, ref_ltag, ref_gene = top_hit.split(':')
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

                        # for logging purposes
                        if feature_is_pseudo:
                            corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']
                            corrected_orf_report[-1][1].qualifiers['pseudo'] = ['']
                        else:
                            corrected_orf_report[-1][0].qualifiers['pseudo'] = ['']

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
                elif 'gene' in feature.qualifiers.keys():
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
                abinit_features_dict.values(),
                seq_ident=seq_ident,
                seq_covg=seq_covg,
                abinit_blast_results=abinit_blast_results,
            )
            prokka_rejects += dropped_abinit_fragments
            logger.debug(f"{seqname}: {len(dropped_abinit_fragments)} gene fragment pairs merged")


            # Check for in-frame conflicts/duplicates
            abinit_features_postprocessed = generate_feature_dictionary(
                abinit_features_postprocessed_list
            )
            logger.info(f"{seqname}: Checking for in-frame overlaps between RATT and ab initio gene annotations")
            unique_abinit_features, inframe_conflicts, abinit_duplicates = find_inframe_overlaps(
                ratt_contig_features,
                abinit_features_postprocessed,
            )
            prokka_rejects += abinit_duplicates
            logger.debug(f"{seqname}: {len(abinit_duplicates)} ab initio ORFs identical to RATT's")
            logger.debug(f"{seqname}: {len(inframe_conflicts)} ab initio ORFs conflicting in-frame with RATT's")
            logger.debug(f'{seqname}: {len(unique_abinit_features)} total ab initio ORFs remain in consideration')


            intergenic_ratt, intergenic_positions, ratt_pre_intergene, ratt_post_intergene = \
                get_interregions(ratt_contig_record_mod, intergene_length=1)
            sorted_intergenic_positions = sorted(intergenic_positions)
            add_features_from_prokka, overlap_conflicts = populate_gaps(
                abinit_features=unique_abinit_features,
                intergenic_positions=sorted_intergenic_positions,
                ratt_pre_intergene=ratt_pre_intergene,
                ratt_post_intergene=ratt_post_intergene,
            )
            add_prokka_contig_record.features += add_features_from_prokka
            logger.debug(f"{seqname}: {len(add_features_from_prokka)} ab initio ORFs fall squarely into RATT's CDS-free regions.")

            # merge the two dicts of lists. hat tip to https://stackoverflow.com/a/5946322
            abinit_conflicts = collections.defaultdict(list)
            for abinit_conflict_group in (inframe_conflicts, overlap_conflicts):
                for key, value in abinit_conflict_group.items():
                    abinit_conflicts[key] += value

            logger.debug(f"{seqname}: {len(abinit_conflicts.keys())} ab initio CDSs in total overlap RATT CDSs. Resolving...")

            for feature_position in abinit_conflicts.keys():
                abinit_feature = unique_abinit_features[feature_position]
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
                    include_abinit, include_ratt, remark = check_inclusion_criteria(
                        ratt_annotation=ratt_feature,
                        abinit_annotation=abinit_feature,
                    )
                    if not include_abinit:
                        prokka_rejects.append((abinit_feature,remark))
                        break
                    # TODO: explain why this is an elif rather than an independent else.
                    elif not include_ratt:
                        ratt_rejects.append((ratt_contig_features_dict.pop(ratt_conflict_loc), remark))
                # Add the abinit feature if it survived all the conflicts
                if include_abinit:
                    add_prokka_contig_record.features.append(abinit_feature)

            for ratt_feature_append in ratt_contig_features_dict.values():
                add_prokka_contig_record.features.append(ratt_feature_append)
            for non_cds in ratt_contig_non_cds:
                add_prokka_contig_record.features.append(non_cds)

            annomerge_records.append(add_prokka_contig_record)

        # Finalize annotation records for this contig
        seqname = '.'.join([isolate_id, contig])
        annomerge_records[i].name = seqname
        # TODO - replace version variable with importlib.version call (and probably url too) in python 3.8+
        annomerge_records[i].annotations['comment'] = "Annotated using hybran " + __version__ + " from https://lpcdrp.gitlab.io/hybran."
        prokka_rec = annomerge_records[i]
        raw_features_unflattened = prokka_rec.features[:]
        raw_features = []
        for f_type in raw_features_unflattened:
            if 'Bio.SeqFeature.SeqFeature' in str(type(f_type)):
                raw_features.append(f_type)
            elif 'list' in str(type(f_type)) and len(f_type) > 0:
                for sub_feature in f_type:
                    if 'Bio.SeqFeature.SeqFeature' not in str(type(sub_feature)):
                        continue
                    else:
                        raw_features.append(sub_feature)
            else:
                continue
        prokka_rec.features = raw_features
        logger.debug(f'{seqname}: final feature annotation verification')
        n_final_cdss = 0
        for feature in prokka_rec.features:
            if feature.type == 'CDS':
                n_final_cdss += 1
                # Adding translated sequences where missing
                if 'translation' not in feature.qualifiers and not designator.is_pseudo(feature.qualifiers):
                    feature_sequence = translate(
                        feature.extract(record_sequence),
                        table=genetic_code,
                        cds=True,
                    )
                    feature.qualifiers['translation'] = [feature_sequence]
                elif designator.is_pseudo(feature.qualifiers):
                    feature.qualifiers.pop('translation', None)
                if 'gene' not in feature.qualifiers.keys():
                    feature.qualifiers['gene'] = feature.qualifiers['locus_tag']

        sorted_final = get_ordered_features(prokka_rec.features)
        prokka_rec.features = sorted_final
        output_isolate_recs.append(prokka_rec)
        isolate_features = prokka_rec.features

        logger.info(f'{seqname}: {n_final_cdss} CDSs annomerge')

    with open(ratt_rejects_logfile, 'w') as ratt_rejects_log:
        [log_feature_fate(_[0], ratt_rejects_log, _[1]) for _ in ratt_rejects]
    with open(prokka_rejects_logfile, 'w') as prokka_rejects_log:
        # we sometimes get the same gene appearing here multiple times, so collapse
        # the list and merge the remarks.
        prokka_rejects.sort(key=lambda _:_[0].qualifiers['locus_tag'][0])
        last_gene = prokka_rejects[0][0]
        last_remark = ''
        for (gene, remark) in prokka_rejects:
            if gene == last_gene:
                if last_remark:
                    remark = ';'.join([last_remark,remark])
            else:
                log_feature_fate(last_gene, prokka_rejects_log, last_remark)
            last_gene = gene
            last_remark = remark
        log_feature_fate(last_gene, prokka_rejects_log, last_remark)

    with open(corrected_abinit_orf_logfile, 'w') as abinit_corlog, \
         open(corrected_ratt_orf_logfile, 'w') as ratt_corlog:
        header = ['locus_tag', 'gene_name', 'strand', 'og_start', 'og_end', 'new_start', 'new_end', 'fixed_start_codon', 'fixed_stop_codon', 'status']
        print('\t'.join(header), file=abinit_corlog)
        print('\t'.join(header), file=ratt_corlog)
        for (orig_feature, corr_feature) in corrected_orf_report:
            if designator.is_raw_ltag(orig_feature.qualifiers['locus_tag'][0]):
                logfile = abinit_corlog
            else:
                logfile = ratt_corlog
            log_coord_correction(orig_feature,
                                 corr_feature,
                                 logfile,
            )

    SeqIO.write(output_isolate_recs, output_genbank, 'genbank')

    logger.debug('annomerge run time: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes')
