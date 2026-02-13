# Title: Merge annotation from Prokka for the positions which are not annotated by RATT
# Description: This program takes as input, a valid EMBL file from RATT annotation or multiple EMBL files in case of
# multiple contigs/chromosome annotations and a Genbank file (.gbk) file from Prokka annotation run with a reference and
# a Genbank file (.gbk) file from Prokka annotation run without a reference. The output is an EMBL file with annotation
# predominantly from RATT and the intergenic regions annotated by RATT are filled with Prokka. This script also
# generates a log file to indicate characteristics of the transferred features from Prokka.

import sys
import collections
import itertools
import os
import tempfile
import pickle
import logging
import time
from copy import deepcopy

# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import networkx as nx

from . import (
    __version__,
    BLAST,
    config,
    converter,
    designator,
    extractor,
)
from .bio import (
    SeqIO,
    sort_features,
    translate,
)
from .compare import (
    compare,
    overlap_inframe,
)
from .demarcate import (
    coord_check,
    has_broken_stop,
)
from .designator import key_ref_gene
from .fusionfisher import fusionfisher
from .lumberjack import (
    log_feature_fates,
    log_coord_corrections,
    log_pseudos,
    log_fusions,
)
from .util import mpbreakpoint

def loc_index(feature):
    """
    Produce the index tuple we use to track feature conflicts.
    """
    return (
        int(feature.location.start),
        int(feature.location.end),
        int(feature.location.strand),
    )

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

    feature.source = extractor.get_source(ref_feature)

def thunderdome(annotation1, annotation2):
    """
    Two genes enter... one gene leaves.
    This function performs a 'standardized' comparison between two
    annotations based on reference-correspondence and pseudo status. Will be used
    in check_inclusion_criteria for the cases dealing with conflicting annotations.
    Note that annotation1 has historically referred to the RATT annotation and annotation2 to the ab initio annotation candidate, so evidence codes/remarks still reference them.

    :param annotation1: SeqFeature
      (annotation1 will be favored if both are found to be equally valid.)
    :param annotation2: SeqFeature
    :returns:
        - include1 (:py:class:`bool`) - whether annotation1 should be kept
        - include2 (:py:class:`bool`) - whether annotation2 should be kept
        - evid (:py:class:`str`) - evidence code for rejected annotation
        - remark (:py:class:`str`) - explanation why the rejected annotation was not included
    """
    logger = logging.getLogger('Thunderdome')
    ann2_delayed_stop = annotation2.de
    ann1_delayed_stop = annotation1.de

    ann2_broken_stop = has_broken_stop(annotation2)[0]
    ann1_broken_stop = has_broken_stop(annotation1)[0]

    if not annotation1.source:
        ann1_coord_status = (True, True)
    else:
        ann1_coord_status = coord_check(annotation1, ref_annotation[
            key_ref_gene(annotation1.source, annotation1.qualifiers['gene'][0])
        ])

    if not annotation2.source:
        ann2_coord_status = (True, True)
    else:
        ann2_coord_status = coord_check(annotation2, ref_annotation[
            key_ref_gene(annotation2.source, annotation2.qualifiers['gene'][0])
        ])

    (ann2_start_ok, ann2_stop_ok) = ann2_coord_status
    ann2_coord_score = sum([int(_) for _ in ann2_coord_status])
    (ann1_start_ok, ann1_stop_ok) = ann1_coord_status
    ann1_coord_score = sum([int(_) for _ in ann1_coord_status])

    ann2_is_pseudo = designator.is_pseudo(annotation2.qualifiers)
    ann1_is_pseudo = designator.is_pseudo(annotation1.qualifiers)

    ann1_longer = (len(annotation2.location) < len(annotation1.location))
    ann2_longer = (len(annotation2.location) > len(annotation1.location))

    if ann2_is_pseudo == ann1_is_pseudo:
        # Both annotations being intact according to their respective reference names
        # suggests that the reference genes are highly similar.
        # RATT's assignment is furthermore based on synteny, so it wins out

        if ((all(ann1_coord_status) and all(ann2_coord_status)) or ann1_coord_status == ann2_coord_status):
            #Both non-pseudo and equal coord_status, use the more complete annotation.
            if not ann2_is_pseudo:
                include2 = False
                include1 = True
                evid = 'forfeit'
                remark = "Equally valid call."
            #Both pseudo
            else:
                if ann2_longer and ann2_delayed_stop and not ann1_delayed_stop:
                    include2 = True
                    include1 = False
                    evid = 'shorter_pseudo'
                    remark = "Rival annotation is longer and contains a valid delayed stop."
                elif ann1_longer and ann1_delayed_stop and not ann2_delayed_stop:
                    include2 = False
                    include1 = True
                    evid = 'shorter_pseudo'
                    remark = "Rival annotation is longer and contains a valid delayed stop."
                #Same location or equivalent delayed stop status
                else:
                    if ann2_broken_stop and not ann1_broken_stop:
                        include2 = False
                        include1 = True
                        evid = 'internal_stop'
                        remark = "Rival annotation doesn't contain any internal stops."
                    elif ann1_broken_stop and not ann2_broken_stop:
                        include2 = True
                        include1 = False
                        evid = 'internal_stop'
                        remark = "Rival annotation doesn't contain any internal stops."
                    else:
                        include2 = False
                        include1 = True
                        evid = 'forfeit'
                        remark = f"Equally valid call."

        elif (all(ann1_coord_status) or ann1_coord_score > ann2_coord_score or (ann1_stop_ok and not ann2_stop_ok)):
            include2 = False
            include1 = True
            evid = 'worse_ref_correspondence'
            remark = "Rival annotation more accurately named and delineated."
        else:
        #This is the only other possibility:
        #elif (all(ann2_coord_status or ann2_coord_score > ann1_coord_score
        #or (ann2_stop_ok and not ann1_stop_ok)):
            include2 = True
            include1 = False
            evid = 'worse_ref_correspondence'
            remark = f"Rival annotation more accurately named and delineated."
    else:
        #Always take the non-pseudo annotation if possible
        if not ann2_is_pseudo and ann1_is_pseudo:
            include2 = True
            include1 = False
            evid = 'pseudo'
            remark = "Non-pseudo annotation takes precedence."
        else:
            include2 = False
            include1 = True
            evid = 'pseudo'
            remark = "Non-pseudo annotation takes precedence."

    if include1 == include2:
        logger.warning(f"Both annotations were marked for {'inclusion' if include1 else 'exclusion'} and one annotation is expected to be excluded:\nDefending Feature\n{annotation1}\n\nChallenging Feature\n{annotation2}\n\nUnhandled scenario."
        )
        include1 = True
        include2 = False
        evid = 'forfeit'
        remark = "Both annotations marked for inclusion. Unhandled scenario."
    return include1, include2, evid, remark


def check_inclusion_criteria(
        annotation1,
        annotation2,
):
    """
    This function compares two annotations and checks for conflicts.
    Either one feature or both will be accepted.
    If there is no conflict, both are kept. Otherwise, they are sent to the thunderdome().

    :param annotation1:
    :param annotation2:
    :returns:
        - include1 (:py:class:`bool`) - whether annotation1 should be kept
        - include2 (:py:class:`bool`) - whether annotation2 should be kept
        - evid (:py:class:`str`) - evidence code for rejected annotation, if any. `None` otherwise.
        - remark (:py:class:`str`) - explanation why the rejected annotation, if any, was not included. `None` otherwise.
    """
    logger = logging.getLogger('CheckInclusionCriteria')
    include1 = True
    include2 = True
    evid = None
    remark = None

    if annotation1.type != annotation2.type:
        pass
    elif annotation1.type != 'CDS':
        # TODO: we should come up with criteria for non-CDS genes
        if annotation1.location == annotation2.location:
            include2 = False
            evid = 'identical_non_cds'
    elif (
            not designator.feature_is_unannotated(annotation1)
            and designator.feature_is_unannotated(annotation2)
    ):
        if overlap_inframe(annotation1.location, annotation2.location):
            include2 = False
            evid = 'unnamed'
            remark = "Unnamed gene and conflicts (overlapping in-frame) with named rival annotation."
    elif (
            not designator.feature_is_unannotated(annotation2)
            and designator.feature_is_unannotated(annotation1)
    ):
        if overlap_inframe(annotation1.location, annotation2.location):
            include1 = False
            evid = 'unnamed'
            remark = "Unnamed gene and conflicts (overlapping in-frame) with named rival annotation."
    else:
        same_gene_name = extractor.get_gene(annotation1) == extractor.get_gene(annotation2)
        same_loc = (annotation1.location == annotation2.location)
        if same_loc or same_gene_name:
            include1, include2, evid, remark = thunderdome(annotation1, annotation2)

        elif not same_gene_name and overlap_inframe(annotation1.location, annotation2.location):
            keepers, fusions, rejects = fusionfisher(
                [annotation1, annotation2],
                ref_annotation,
                adjudicate=False,
            )
            for trial in rejects:
                reject = trial['feature']
                if reject == annotation1:
                    include1 = False
                    evid = trial['evid']
                    remark = trial['remark']
                elif reject == annotation2:
                    include2 = False
                    evid = trial['evid']
                    remark = trial['remark']
            # fusionfisher didn't detect a misannotation, but it didn't detect a fusion either.
            # welcome to the thunderdome!
            if not fusions and not rejects:
                include1, include2, evid, remark = thunderdome(annotation1, annotation2)

        else:
            #include everything if different names and not overlapping in frame
            include1 = True
            include2 = True

    return include1, include2, evid, remark

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

def merge(overlap_G):
    """
    Given an overlap graph, check for conflicting annotations and choose the "best" ones, returning a full set of consistent features.
    :param overlap_G:
      a networkx Graph as produced by compare() or cross_examine().
      Nodes in this graph must have an 'annotation' attribute containing the corresponding SeqFeature object.
      The SeqFeature object should have a label attribute corresponding to the node label.
    :return: list of sorted SeqFeatures being the result of the merge
    """
    refined_G = overlap_G.copy()
    rejects_data = []

    for cc in nx.connected_components(overlap_G):
        cc_features = [overlap_G.nodes[n]['annotation'] for n in cc]

        # singleton connected components are standalone/unique features
        if len(cc_features) == 1:
            continue

        for f1, f2 in itertools.combinations(cc_features, 2):
            # skip comparison if one of the contenders has been previously defeated
            if f1.label not in refined_G or f2.label not in refined_G:
                continue

            include_f1, include_f2, evidence, remark = check_inclusion_criteria(f1, f2)
            if not include_f1:
                refined_G.remove_node(f1.label)
                rejects_data.append({
                            'feature': f1,
                            'superior': f2,
                            'evid': evidence,
                            'remark':remark,
                })
            if not include_f2:
                refined_G.remove_node(f2.label)
                rejects_data.append({
                            'feature': f2,
                            'superior': f1,
                            'evid': evidence,
                            'remark':remark,
                })

    return (
        sort_features([refined_G.nodes[n]['annotation'] for n in refined_G]),
        rejects_data,
    )

def run(
        isolate_id,
        organism,
        strain,
        genome,
        annotation_fp,
        ref_proteins_fasta,
        ref_gbk_list,
        script_directory,
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
    ref_annotation = extractor.load_gbks(ref_gbk_list, feature_types=['CDS'])

    annomerge_records = list(SeqIO.parse(genome, "fasta"))
    contigs = [record.id for record in annomerge_records]

    ratt_file_path = os.path.join(file_path, 'ratt')
    ratt_features = ratt.postprocess(
        isolate_id,
        contigs,
        ratt_outdir=ratt_file_path,
        postprocess_outdir=os.path.join(file_path, 'ratt-postprocessed'),
        ref_annotation=ref_annotation,
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
            ) = compare(ratt_contig_features, abinit_contig_features, eliminate_colocated=False)

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
                    overlap_G.remove_node(abinit_feature.label)
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
                # Skip feature if it was previously rejected during check of co-located features.
                #     Setting compare()'s eliminate_colocated=True would avoid this,
                #     but it would also exclude those colocated RATT features from being considered
                #     as conflicting with anything else.
                #     For identifying fusions, this isn't desirable.
                if abinit_feature.label not in overlap_G.nodes:
                    continue

                for _, ratt_feature_label in overlap_G.edges(abinit_feature.label):
                    ratt_feature = overlap_G.nodes[ratt_feature_label]['annotation']
                    abinit_conflicts[loc_index(abinit_feature)].append(loc_index(ratt_feature))

            logger.info(f"{seqname}: {len(inframe_conflicts)} ab initio ORFs conflicting in-frame with RATT's")
            logger.info(f"{seqname}: {len(abinit_conflicts)} ab initio features in total overlap RATT features. Resolving...")


            for feature_position in abinit_conflicts.keys():
                abinit_feature = abinit_contig_features_dict[feature_position]
                # Conflict Resolution
                for ratt_conflict_loc in abinit_conflicts[feature_position]:
                    # if the RATT annotation got rejected at some point, its remaining conflicts are moot
                    if ratt_conflict_loc not in ratt_contig_features_dict.keys():
                        include_abinit = True
                        continue
                    ratt_feature = ratt_contig_features_dict[ratt_conflict_loc]
                    include_ratt, include_abinit, evid, remark = check_inclusion_criteria(
                        ratt_feature,
                        abinit_feature,
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
                    feature.qualifiers['gene'] = [ feature.qualifiers['locus_tag'][0] ]

        sorted_final = sort_features(annomerge_records[i].features)
        annomerge_records[i].features = sorted_final

        logger.info(f'{seqname}: {n_final_cdss} CDSs annomerge')

    with open(ratt_rejects_logfile, 'w') as ratt_rejects_log:
        log_feature_fates(ratt_rejects, logfile=ratt_rejects_log)
    with open(prokka_rejects_logfile, 'w') as prokka_rejects_log:
        log_feature_fates(prokka_rejects, logfile=prokka_rejects_log)

    annomerge_records_dict = {i: annomerge_records[i].features for i in range(len(annomerge_records))}
    corrected_orf_logfile = os.path.join(
        isolate_id,
        'annomerge',
        'coord_corrections.tsv'
    )
    with open(corrected_orf_logfile, 'w') as corr_log:
        log_coord_corrections(annomerge_records_dict, corr_log)


    designator.assign_locus_tags(
        annomerge_records_dict,
        prefix=isolate_id,
    )
    SeqIO.write(annomerge_records, output_genbank, "genbank")

    pseudoscan_logfile = os.path.join(
        isolate_id,
        'pseudoscan_report.tsv'
    )
    with open(pseudoscan_logfile, 'w') as p_log:
        log_pseudos(annomerge_records_dict, p_log)

    fusionfisher_logfile = os.path.join(
        isolate_id,
        'fusion_report.tsv'
    )
    with open(fusionfisher_logfile, 'w') as f_log:
        log_fusions(annomerge_records_dict, f_log)

    logger.debug('postprocessing and annomerge run time: ' + str(int((time.time() - start_time) / 60.0)) + ' minutes')
