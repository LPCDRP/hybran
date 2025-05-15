from copy import deepcopy
import logging

from . import (
    designator,
    extractor,
    pseudoscan,
)
from .annomerge import merge_qualifiers
from .bio import FeatureProperties
from .compare import (
    complementary,
    overlap_inframe,
)
from .demarcate import (
    coord_check,
    update_termini,
)
from .designator import key_ref_gene


def fissionfuser(flist, ref_annotation, seq_ident, seq_covg):
    """
    Given a list of features ordered by genomic position, identify adjacent gene fragments and combine them into a single feature.
    :param flist: list of SeqFeature objects
    :param ref_annotation: dict mapping the reference IDs and gene names (keys according to annomerge.key_ref_gene) to (Autarkic)SeqFeatures
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
