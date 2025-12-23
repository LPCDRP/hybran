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
from copy import deepcopy

# standard multiprocessing can't pickle lambda
import multiprocess as multiprocessing
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

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
            keepers, fusions, rejects = fusionfisher(
                [ratt_annotation, abinit_annotation],
                ref_annotation,
                adjudicate=False,
            )
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
