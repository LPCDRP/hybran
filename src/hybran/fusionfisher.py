from copy import deepcopy
import logging

from . import (
    designator,
    extractor,
    pseudoscan,
)
from .compare import (
    complementary,
    have_same_stop,
)
from .demarcate import (
    coord_check,
    update_termini,
)
from .designator import key_ref_gene


def fusion_name(feature1, feature2):
    """
    For use by fusionfisher.
    Create a fusion gene and product name. Usually, it's just "gene1::gene2" and "product1 / product2",
    but we need to watch out for redundancies, especially when one of the two features is already a fusion
    including the second's name.
    """
    f1_name_components = extractor.get_gene(feature1).split('::')
    f2_name_components = extractor.get_gene(feature2).split('::')
    reduced_f2_name_components = [
        comp_name for comp_name in f2_name_components
        if comp_name not in f1_name_components
    ]
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
    reduced_f2_product_components = [
        comp_prod for comp_prod in f2_product_components
        if comp_prod not in f1_product_components
    ]
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
        base.qualifiers['product'] = [product]

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

def add_components(orig_components, new_components):
    """
    Add fusion components non-redundantly
    :param orig_component: list of SeqFeatures in the existing component list
    :param new_components: list of SeqFeatures to add
    :return: list of SeqFeatures representing all components non-redudantly
    """

    if not orig_components:
        final_components = new_components
    else:
        orig_component_gene_names = [extractor.get_gene(f) for f in orig_components]
        new_component_gene_names  = [extractor.get_gene(f) for f in new_components]

        final_components = orig_components.copy()
        for i in range(len(new_component_gene_names)):
            if (
                    # Do not add components if they represent the same gene as an existing component
                    #   (we're not handling the case of multiple copies of a gene getting fused together with something else...
                    new_component_gene_names[i] in orig_component_gene_names
                    # And don't list a fusion gene itself as a fusion component
                    or new_components[i].fusion_type
            ):
                continue
            final_components.append(new_components[i])
        final_components = sorted(final_components, key=lambda f: (f.location.start, f.location.end))

    return final_components

def fusionfisher(feature_list, ref_annotation, adjudicate=True):
    """
    This function parses through a list of CDSs and returns a unique list of CDSs, cleaning up annotation artifacts due to gene fusion events, as well as renaming such genes and those that have conjoined with their neighbor.

    :param feature_list: list of sorted SeqFeature objects.
    :param ref_annotation: dict mapping the reference IDs and gene names (keys according to annomerge.key_ref_gene) to (Autarkic)SeqFeatures
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
        'partial': [],
        'whole': [],
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
            # Artifact due to a partial gene fusion
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

                component1 = deepcopy(upstream)
                component2 = deepcopy(downstream)
                coord_check(
                    component1,
                    ref_annotation[key_ref_gene(component1.source, extractor.get_gene(component1))],
                    fix_stop=True,
                    seek_stop=False,
                    check_context=False,
                    best_effort=True,
                )
                coord_check(
                    component2,
                    ref_annotation[key_ref_gene(component2.source, extractor.get_gene(component2))],
                    fix_start=True,
                    seek_stop=False,
                    check_context=False,
                    best_effort=True,
                )

                fusion_upgrade(
                    base=prev_feature,
                    upstream=upstream,
                    downstream=downstream,
                )

                rejects.append({
                    'feature':outlist.pop(),
                    'superior':prev_feature,
                    'evid':'combined_annotation',
                    'remark':"Apparent partial fusion gene. Name incorporated into rival feature's and redundant locus removed.",
                })
                remarkable['partial'].append(prev_feature)
                if not prev_feature.fusion_type:
                    prev_feature.fusion_type = 'partial'
                prev_feature.fusion_components = add_components(
                    prev_feature.fusion_components,
                    [component1, component2],
                )

        elif have_same_stop(prev_feature.location, feature.location):
            #
            # whole gene fusions
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

                component1 = deepcopy(upstream)
                coord_check(
                    component1,
                    ref_annotation[key_ref_gene(component1.source, extractor.get_gene(component1))],
                    fix_stop=True,
                    seek_stop=False,
                )
                fusion_upgrade(
                    base=upstream,
                    upstream=upstream,
                    downstream=downstream,
                )

                remarkable['whole'].append(upstream)
                # A redundant partial gene fusion would look like a whole gene fusion
                # see 'idempotence' test case.
                if not upstream.fusion_type:
                    upstream.fusion_type = 'whole'
                upstream.fusion_components = add_components(
                    upstream.fusion_components,
                    [component1, downstream],
                )
            #
            # Another signature of a partial gene fusion
            #
            elif complementary(prev_feature, feature):
                if feature.rcs:
                    upstream, downstream = feature, prev_feature
                else:
                    upstream, downstream = prev_feature, feature

                component1 = deepcopy(upstream)
                component2 = deepcopy(downstream)
                coord_check(
                    component1,
                    ref_annotation[key_ref_gene(component1.source, extractor.get_gene(component1))],
                    fix_stop=True,
                    seek_stop=False,
                    check_context=False,
                    best_effort=True,
                )
                coord_check(
                    component2,
                    ref_annotation[key_ref_gene(component2.source, extractor.get_gene(component2))],
                    fix_start=True,
                    seek_stop=False,
                    check_context=False,
                    best_effort=True,
                )
                fusion_upgrade(
                    base=upstream,
                    upstream=upstream,
                    downstream=downstream,
                    update_location=True,
                )

                remarkable['partial'].append(upstream)
                if not upstream.fusion_type:
                    upstream.fusion_type = 'partial'
                upstream.fusion_components = add_components(
                    upstream.fusion_components,
                    [component1, component2],
                )
                rejects.append({
                    'feature':downstream,
                    'superior':upstream,
                    'evid':'combined_annotation',
                    'remark':"Apparent partial gene fusion.",
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


    return outlist, remarkable['partial'] + remarkable['whole'], rejects
