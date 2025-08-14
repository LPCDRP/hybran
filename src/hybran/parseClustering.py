import os
import re
import tempfile
import subprocess
import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

from . import (
    BLAST,
    CDHIT,
    config,
    designator,
    extractor,
)
from .bio import SeqIO


def parse_clustered_proteins(clustered_proteins, annotations):
    """
    Parses clustered_proteins and categorizes each cluster into four categories
    1. Clusters that have unannotated gene names and annotated gene names that are the same
    2. Clusters that have a single gene
    3. Clusters that have unannotated gene names and annotated gene names that are difference
    4. Clusters that only have unannotated gene names

    :param clustered_proteins: file name of clustering output ('clustered_proteins')
    :param annotations: directory that has the annotations used to create the clusters
    :return: list of dictionaries
    """

    # Private Function #
    ##################################################################################################################
    def gff_dict(gffs):
        """
        Parses all GFFs in the annotation_dir and gets gene names and locus tags.
        Creates a ditionary of dictionaries with isolate ID as the first key
        and GFF ID as the inner key. Locus tag and gene names are tuple values
        for the inner key

        :param gffs: str directory of GFFs
        :return: dictionary of dictionaries
        """
        gff_dictionary = {}
        for gff in gffs:
            isolate_id_ltag = {}
            try:
                isolate_id = os.path.splitext(os.path.basename(gff))[0]
                if gff.endswith('.gff'):
                    with open(gff, 'r') as gff_file:
                        for line in gff_file:
                            if line.startswith('#'):
                                continue
                            gene = ''
                            column = line.rstrip('\n').split('\t')
                            if len(column) >= 8:
                                info = column[8].split(';')
                                gene_id = ''.join([i.split('=')[1] for i in info if i.startswith('ID=')])
                                locus_tag = ''.join([i.split('=')[1] for i in info if i.startswith('locus_tag=')])
                                if not gene_id or not locus_tag:
                                    continue
                                gene = ','.join([i.split('=')[1] for i in info if i.startswith('gene=')])
                                isolate_id_ltag[gene_id] = (locus_tag, gene)
                        if not isolate_id in gff_dictionary.keys():
                            gff_dictionary[isolate_id] = isolate_id_ltag
                        # this is the case if the reference genome itself is being processed as a sample,
                        # which one would do if wanting the annotation updated with ab initio predictions in the gaps
                        else:
                            gff_dictionary[isolate_id].update(isolate_id_ltag)
            except IOError:
                continue
        return gff_dictionary

    # Different gene names/references in cluster
    def different_genes(gene_cluster_list):
        """
        Collates all unique gene names and locus tags in a given cluster

        :param gene_cluster_list: list of tuples containing locus tags and gene names
        :return: two lists
        """
        genes = list(set(sorted([gene[1] for gene in gene_cluster_list])))
        locus_tags = list(set(sorted([locus[0] for locus in gene_cluster_list if locus[0]])))
        return genes, locus_tags

    ##################################################################################################################

    gffs = gff_dict(annotations)
    representative_fasta_list = []
    gene_cluster = OrderedDict()
    different_genes_cluster_w_ltags = {}
    same_genes_cluster_w_ltags = {}
    underscores = {}
    l_tag_only_clusters = {}
    unique_genes_list = []
    with open(clustered_proteins, 'r') as clustered_proteins_file:
        for line in clustered_proteins_file:
            cluster_members = [_.groups() for _ in re.finditer(r'(?P<sample_id>\S+)@@@(?P<seq_id>\S+)', line)]
            (rep_isolate, rep_seq_id) = cluster_members[0]
            representative_ltag_gene_tup = gffs[rep_isolate][rep_seq_id]
            representative = ','.join([rep_isolate, ','.join(representative_ltag_gene_tup)])

            representative_fasta_list.append([rep_isolate] + list(representative_ltag_gene_tup))
            gene_cluster[representative] = []

            cluster_list = []
            cluster_list_w_isolate = []
            cluster_list.append(representative_ltag_gene_tup)
            cluster_list_w_isolate.append([rep_isolate] + list(representative_ltag_gene_tup))
            for isolate, gene_id in cluster_members[1:]:
                cluster_list.append(gffs[isolate][gene_id])
                cluster_list_w_isolate.append([isolate] + list(gffs[isolate][gene_id]))
                gene_cluster[representative].append(','.join([isolate, ','.join(gffs[isolate][gene_id])]))
            number_genes = list(set(sorted([gene for locus_tag, gene in cluster_list])))
            unique_genes = []
            for data in number_genes:
                if not designator.is_raw_ltag(data):
                    unique_genes.append(data)
            # Getting clusters with only 1 gene that is an L tag or an underscore
            if (
                    len(unique_genes) == 1
                    and any([designator.is_raw_ltag(gene) for locus_tag, gene in cluster_list])
            ):
                same_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            elif (len(cluster_members) == 1 and
                    (designator.is_raw_ltag(representative_ltag_gene_tup[1]))):
                unique_genes_list.append([rep_isolate] + list(representative_ltag_gene_tup))
            # All different with L tags
            elif (
                    len(unique_genes) > 1
                    and any([designator.is_raw_ltag(gene) for locus_tag, gene in cluster_list])
            ):
                different_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            #
            # L tag only clusters
            elif all(designator.is_raw_ltag(gene) for locus_tag, gene in cluster_list):
                l_tag_only_clusters[representative] = cluster_list_w_isolate
    return [different_genes_cluster_w_ltags, same_genes_cluster_w_ltags, l_tag_only_clusters,
            unique_genes_list]


def get_gene_name(id_string):
    if '|' in id_string:
        [seqid, locus_tag, gene_name] = id_string.split('|')
        if gene_name:
            return gene_name
        else:
            return locus_tag
    else:
        return id_string

def check_matches_to_known_genes(
        query_seq,
        reference_seqs,
        generic_seqs,
        cluster_type,
        orf_increment,
        seq_ident,
        seq_covg,
):
    """
    Use BLAST to search for any query_seq hits to reference_seqs or,
    failing that, generic_seqs meeting the seq_ident and seq_covg thresholds.
    query_seq's id and description attributes also get updated and its record
    appended to the generic_seqs fasta file if no adequate match to known
    genes is found.

    :param query_seq: SeqRecord query gene AA sequence
    :param reference_seqs: str fasta file name of reference genes to map against
    :param generic_seqs: str fasta file name of hitherto assigned generic genes (ORF####)
    :param cluster_type: str category of MCL clusters that the query_seq came from
    :param orf_increment: int increment to start new ORFs
    :returns:
        - name_to_assign (:py:class:`str`)
        - query_gene_closest_refs (:py:class:`list`) -
              list of tab delimited strings detailing closest matches
              among reference and generic seqs. Empty list if a
              match was successfully found and a new generic name
              was not assigned
        - orf_increment (:py:class:`int`) - updated increment of ORFs
    """
    best_subcriticals = []
    query_gene_closest_refs = []
    assign_new_generic = True

    if os.stat(generic_seqs).st_size:
        refs = [reference_seqs, generic_seqs]
    else:
        refs = [reference_seqs]

    for ref in refs:
        if cluster_type == 'multiref':
            top_hit, blast_stats = BLAST.bidirectional_best_hit(
                query=query_seq,
                subject=ref,
                min_bitscore=config.cnf.blast.min_bitscore,
                min_seq_covg=config.cnf.blast.min_coverage,
                identify=get_gene_name,
                #strict=True,
            )
            if top_hit:
                name_to_assign = top_hit[query_seq.id]
                assign_new_generic = False
        else:
            top_hit = None

        #if not top_hit:
        #    best_subcriticals += check_subcriticals(blast_stats)
        best_subcriticals = []

    if assign_new_generic:
        (orf_id, orf_increment) = designator.assign_orf_id(orf_increment)
        query_seq.id = orf_id
        query_seq.description = 'False'
        with open(generic_seqs, 'a') as orf_fasta:
            SeqIO.write(query_seq, orf_fasta, 'fasta')
        name_to_assign = orf_id
        # prepare the novelty report data.
        #
        # if we didn't get *any* hits at all, we'll have a basically null row
        # here.
        if not best_subcriticals:
            query_gene_closest_refs = ['\t'.join([cluster_type, orf_id])]
        else:
            query_gene_closest_refs = ['\t'.join([cluster_type, orf_id, hit])
                                       for hit in best_subcriticals]

    return name_to_assign, query_gene_closest_refs, orf_increment

def check_subcriticals(subcriticals):
    """
    Check hits that didn't pass thresholds and find the best among them

    :param subcriticals: dict from BLAST.summarize() for remaining hits
    :returns:
       - empty list (if there are no hits at all) or list of three report strings
         for the best subcritical matches by each category.
    """
    if not subcriticals:
        return []
    best_subc_iden = BLAST.top_hit(subcriticals, metric='iden')
    best_subc_scov = BLAST.top_hit(subcriticals, metric='scov')
    best_subc_qcov = BLAST.top_hit(subcriticals, metric='qcov')
    best_subcriticals = [
        '\t'.join([_[0],
                   _[1],
                   str(subcriticals[_[0]]['iden']),
                   str(subcriticals[_[0]]['scov']),
                   str(subcriticals[_[0]]['qcov'])])
        for _ in zip([best_subc_iden, best_subc_scov, best_subc_qcov],
                     ['identity','scov','qcov'])
    ]
    return best_subcriticals


def update_dictionary_ltag_assignments(isolate_id, isolate_ltag, new_gene_name):
    """
    Updates a global dictionary with new gene name assignments

    :param isolate_id: str isolate ID
    :param isolate_ltag: str locus tag to update
    :param new_gene_name: str gene name to assign
    :return: None
    """
    logger = logging.getLogger('NameGene')
    # TODO: not doing pseudogene calling on name assignments from clustering.
    #       Instead, we're requiring full-coverage matches at this point until
    #       we can work out a way to take care of the fallout of a pseudo name
    #       assignment (checking neighboring genes, potentially reevaluating conflicting annotations...)
    pseudo = False

    if isolate_id not in isolate_update_dictionary.keys():
        isolate_update_dictionary[isolate_id] = {}
        isolate_update_dictionary[isolate_id][isolate_ltag] = dict(
            name = new_gene_name,
            pseudo = pseudo,
        )
        logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
    else:
        isolate_dict_added = isolate_update_dictionary[isolate_id]
        if isolate_ltag not in isolate_dict_added.keys():
            logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
            isolate_update_dictionary[isolate_id][isolate_ltag] = dict(
                name = new_gene_name,
                pseudo = pseudo,
            )
    return


def get_cluster_fasta(rep, cluster_list):
    """
    Identifies all remaining unannotated genes and adds them to
    the unannotated sequence FASTA

    :param rep: str gene name
    :param cluster_list: list of isolates that have the gene
    :return: FASTA file name and list of unannotated genes
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    rep_isolate_id = rep.split(',')[0]
    rep_locus = rep.split(',')[1]
    rep_gene_name = rep.split(',')[2]
    cluster_fasta_fp = tempfile.NamedTemporaryFile(suffix='.fasta',
                                                   dir=hybran_tmp_dir,
                                                   delete=False, mode='w')
    cluster_temp_fasta = cluster_fasta_fp.name
    added_seq = []
    unannotated_genes_list = []
    fasta_tmp = open(cluster_temp_fasta, 'w')
    if not designator.is_raw_ltag(rep_gene_name):
        rep_sequence = isolate_sequences[rep_isolate_id][rep_locus]
        rep_sequence_object = rep_sequence
        seq_id = '|'.join([rep_isolate_id, rep_locus, rep_gene_name])
        added_seq.append(seq_id)
        rep_record = SeqRecord(rep_sequence_object, id=seq_id)
        SeqIO.write(rep_record, fasta_tmp, 'fasta')
    else:
        rep_list = [rep_isolate_id, rep_locus, rep_gene_name]
        unannotated_genes_list.append(rep_list)
    for gene in cluster_list:
        gene_isolate_id, gene_locus, gene_name = gene
        if designator.is_raw_ltag(gene_name):
            unannotated_genes_list.append(gene)
        else:
            gene_id = '|'.join([gene_isolate_id, gene_locus, gene_name])
            if gene_id in added_seq:
                continue
            else:
                gene_sequence = isolate_sequences[gene_isolate_id][gene_locus]
                gene_sequence_object = gene_sequence
                added_seq.append(gene_id)
                gene_record = SeqRecord(gene_sequence_object, id=gene_id)
                SeqIO.write(gene_record, fasta_tmp, 'fasta')
    return cluster_temp_fasta, unannotated_genes_list


def cluster_annotation_presence(cluster_list):
    """
    Checks if all genes in the cluster list have a gene name

    Assumes representative sequence is also present in cluster_list
    :param cluster_list: List of tuples
    :return: bool
    """

    annotated = 0
    for gene in cluster_list:
        gene_locus = gene[1]
        gene_name = gene[2]
        if not designator.is_raw_ltag(gene_name):
            annotated += 1
    return (annotated == len(cluster_list))


def unique_seqs(annotations):
    """
    Creates a FASTA with unique sequences based on a given directory
    of reference Genbank files. Uses GFFs from referene and input
    genomes to create content for the output ref_cdss_protein-all.fasta

    :param annotations: list of annotation files in GFF format
    :return: str FASTA file name and a dictionary of dictionaries
    """
    hybran_tmp_dir = config.hybran_tmp_dir

    def grep_seqs(gff):
        """
        Finds all translations in a given GFF

        :param gff: str GFF file anme
        :return: list of GFF lines
        """
        gff_lines = []
        cmd = ['grep', 'translation=', gff]
        translations = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        for line in translations.stdout:
            try:
                gff_lines.append(line)
            except TypeError:
                continue
        return gff_lines
    protein_cds = []
    isolate_seqs = {}
    for gff in annotations:
        if gff.endswith('.gff'):
            raw_out = grep_seqs(gff)
            gff_name = os.path.splitext(os.path.basename(gff))[0]
            if not gff_name in isolate_seqs.keys():
                isolate_seqs[gff_name] = {}
            for line in raw_out:
                if line.split('\t')[2] == 'CDS':
                    gene = None
                    for i in line.split(';'):
                        if i.startswith('gene='):
                            gene = [i.split('=')[1].rstrip('\n')][0]
                        if i.startswith('locus_tag='):
                            locus_tag = [i.split('=')[1].rstrip('\n')][0]
                    if not gene:
                        gene = locus_tag
                    translation = [i.split('=')[1] for i in line.split(';') if i.startswith('translation=')][0]
                    eggnog_annot = [i for i in line.split(';') if i.startswith('note=') and 'Eggnog' in i.split('=')[1]]
                    eggnog = False
                    if eggnog_annot:
                        eggnog = True
                    record = SeqRecord(Seq(translation.rstrip('\n')),
                                       id=gene,
                                       description=str(eggnog) + '|' + gff_name)
                    isolate_seqs[gff_name][locus_tag] = Seq(translation.rstrip('\n'))
                    protein_cds.append(record)
    all_proteins = hybran_tmp_dir + '/clustering/cdss_protein-all.fasta'
    final_proteins = hybran_tmp_dir + '/clustering/cdss_protein.fasta'
    # Write ref_cdss_protein-all.fasta
    with open(all_proteins, 'w') as ref_cds_out:
        for s in protein_cds:
            SeqIO.write(s, ref_cds_out, 'fasta')
    # Run CD-HIT with 99% thresholds
    CDHIT.run(nproc=1,
              fasta=all_proteins,
              seq_ident=99, seq_covg=99,
              out=final_proteins)
    return final_proteins, isolate_seqs

def prepare_for_eggnog(unannotated_seqs, outfile):
    """
    Creates a FASTA file of unannotated sequences that have not
    been annotated through EggNOG

    :param unannotated_seqs: str FASTA file name
    :param outfile: file handle for output FASTA file name
    :return: None
    """
    eggnog_seqs = []
    for record in SeqIO.parse(unannotated_seqs, 'fasta'):
        header = record.description.split('|')
        state = header[0]
        if 'False' in state:
            eggnog_seqs.append(record)
    SeqIO.write(eggnog_seqs, outfile, 'fasta')


def singleton_clusters(singleton_dict, reference_fasta, unannotated_fasta, orf_increment, seq_ident, seq_covg):
    """
    If a cluster has only one gene (with no gene names),
    then BLAST representative sequence to the reference and if there
    is a hit with specified amino acid and coverage thresholds,
    all candidate novel genes in the cluster is annotated
    with the reference gene. If the representative does not hit a reference,
    assign a generic name to the genes in the  cluster.

    :param singleton_dict: dict of clusters
    :param reference_fasta: str reference FASTA proteome
    :param unannotated_fasta: str unannotated FASTA file
    :param orf_increment: int increment to start new ORFs
    :return: int updated increment of ORFs
    """
    logger = logging.getLogger('SingletonClusters')

    logger.debug('Number of singleton clusters with single genes: ' + str(len(singleton_dict)))

    out_list = []
    novel_genes_closest_refs = []
    for single_gene in singleton_dict:
        isolate_id = single_gene[0]
        locus_tag = single_gene[1]
        gene_name = single_gene[2]
        out_list.append(str(single_gene))
        if designator.is_raw_ltag(gene_name):
            gene_sequence = isolate_sequences[isolate_id][locus_tag]
            name_to_assign, query_closest_refs, orf_increment = check_matches_to_known_genes(
                query_seq = SeqRecord(gene_sequence, id=gene_name),
                reference_seqs = reference_fasta,
                generic_seqs = unannotated_fasta,
                cluster_type = 'singleton',
                orf_increment = orf_increment,
                seq_ident = seq_ident,
                seq_covg = seq_covg,
            )
            novel_genes_closest_refs += query_closest_refs
            update_dictionary_ltag_assignments(isolate_id, locus_tag, name_to_assign)
        else:
            continue
    with open('./clustering/singleton_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return orf_increment, novel_genes_closest_refs


def only_ltag_clusters(in_dict, reference_fasta, unannotated_fasta, orf_increment, seq_ident, seq_covg):
    """
    If a cluster has only candidate novel genes (with no gene names),
    then BLAST representative sequence to reference
    and if there is a hit with specified amino acid
    and coverage thresholds, all candidate novel genes in the cluster
    is annotated with the reference gene. If the representative does
    not hit a reference gene, assign a ORF locus tag to the genes
    in the cluster.

    :param in_dict: dict of clusters
    :param reference_fasta: str reference FASTA proteome
    :param unannotated_fasta: str unannotated FASTA file
    :param orf_increment: int increment to start new ORFs
    :return: int updated increment of ORFs
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('NewGeneClusters')

    logger.debug('Number of clusters with no gene names: ' + str(len(list(in_dict.keys()))))

    out_list = []
    novel_genes_closest_refs = []
    for rep_gene in in_dict.keys():
        out_list.append(rep_gene + ': ' + str(in_dict[rep_gene]))
        isolate = rep_gene.split(',')[0]
        locus = rep_gene.split(',')[1]
        gene = rep_gene.split(',')[2]
        rep_sequence = isolate_sequences[isolate][locus]
        rep_record = SeqRecord(rep_sequence, id='', description='')
        name_to_assign, query_closest_refs, orf_increment = check_matches_to_known_genes(
            query_seq = rep_record,
            reference_seqs = reference_fasta,
            generic_seqs = unannotated_fasta,
            cluster_type = 'no_refs',
            orf_increment = orf_increment,
            seq_ident = seq_ident,
            seq_covg = seq_covg,
        )
        novel_genes_closest_refs += query_closest_refs
        update_dictionary_ltag_assignments(isolate, locus, name_to_assign)
        for other_genes_in_cluster in in_dict[rep_gene]:
            update_dictionary_ltag_assignments(
                other_genes_in_cluster[0],
                other_genes_in_cluster[1],
                name_to_assign,
            )
    with open('./clustering/onlyltag_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return orf_increment, novel_genes_closest_refs


def single_gene_clusters(single_gene_dict):
    """
    Handling clusters with only one reference/gene name in cluster
    If cluster has candidate novel genes (L_***** or L2_*****)
    clustered together with a reference annotation, update
    all candidate novel genes in this cluster with the reference gene name.

    :param single_gene_dict: dict of clusters
    :return: None
    """
    logger = logging.getLogger('ClustersWithOneGeneName')

    rep_ltag_keys = []
    logger.debug('Number of clusters with single genes: ' + str(len(list(single_gene_dict.keys()))))

    for key_sgc in single_gene_dict:
        rep_ltag = False
        (rep_sample, rep_locus_tag, rep_gene_name) = key_sgc.split(',')
        if designator.is_raw_ltag(rep_gene_name):
            rep_ltag = True
            rep_ltag_keys.append(key_sgc)
        else:
            gene_to_add = rep_gene_name
        # If representative sequence is a reference gene and not an L-tag
        if not rep_ltag:
            genes_in_cluster = single_gene_dict[key_sgc]
            for sample, locus_tag, gene_name in genes_in_cluster:
                if gene_name and designator.is_raw_ltag(gene_name):
                    update_dictionary_ltag_assignments(sample, locus_tag, gene_to_add)
        # If representative sequence is a L-tag
        else:
            genes_in_cluster = single_gene_dict[key_sgc]
            gene_to_add = ''
            for sample, locus_tag, gene_name in genes_in_cluster:
                if gene_name and not designator.is_raw_ltag(gene_name):
                    gene_to_add = gene_name
                    break
            if gene_to_add:
                for sample, locus_tag, gene_name in genes_in_cluster:
                    if designator.is_raw_ltag(gene_name):
                        update_dictionary_ltag_assignments(sample, locus_tag, gene_to_add)


def multigene_clusters(in_dict, single_gene_cluster_complete, unannotated_fasta, orf_increment, seq_ident, seq_covg):
    """
    If a cluster has multiple reference genes and L_tags,
    BLAST the L_tags to all genes in the cluster that are reference
    and annotate L_tag with the top hit.

    :param in_dict: dict of clusters
    :param single_gene_cluster_complete: dict of clusters that have a single sequence
    :param unannotated_fasta: str unannotated FASTA file
    :param orf_increment: int increment to start new ORFs
    :return: int updated increment of ORFs
             list report of closest matches for novel genes
    """
    logger = logging.getLogger('ClustersWithManyGeneNames')
    logger.debug('Number of clusters that have many gene names and genes with no name: ' + str(len(list(in_dict.keys()))))

    num_multi = 0
    out_list = []
    novel_genes_closest_refs = []
    for gene in in_dict.keys():
        num_multi += 1
        unassigned_l_tags = []
        true_multi_cluster = False
        gene_elements = gene.split(',')
        out_list.append(gene + ': ' + str(in_dict[gene]))
        for gene_in_cluster in in_dict[gene]:
            (sample, locus_tag, gene_name) = gene_in_cluster
            if gene_name in gene_elements:
                continue
            elif designator.is_raw_ltag(gene_name):
                unassigned_l_tags.append(gene_in_cluster)
            else:
                true_multi_cluster = True
        if true_multi_cluster:
            all_genes_annotated = cluster_annotation_presence(in_dict[gene])
            # If all genes in the cluster are annotated, do nothing. If there is an unannotated gene in the cluster,
            # blast it to all annotated genes in cluster and annotated with top hit
            if not all_genes_annotated:
                fasta_fp_to_blast, genes_to_annotate = get_cluster_fasta(gene, in_dict[gene])
                for unannotated_gene in genes_to_annotate:
                    unannotated_gene_isolate = unannotated_gene[0]
                    unannotated_gene_locus = unannotated_gene[1]
                    unannotated_gene_seq = SeqRecord(isolate_sequences[unannotated_gene_isolate][unannotated_gene_locus], id=unannotated_gene_locus)
                    name_to_assign, query_closest_refs, orf_increment = check_matches_to_known_genes(
                        query_seq = unannotated_gene_seq,
                        reference_seqs = fasta_fp_to_blast,
                        generic_seqs = unannotated_fasta,
                        cluster_type = 'multiref',
                        orf_increment = orf_increment,
                        seq_ident = seq_ident,
                        seq_covg = seq_covg,
                    )
                    novel_genes_closest_refs += query_closest_refs
                    update_dictionary_ltag_assignments(
                        unannotated_gene_isolate,
                        unannotated_gene_locus,
                        name_to_assign,
                    )
        else:
            if len(unassigned_l_tags) > 0:
                single_gene_cluster_complete[gene] = in_dict[gene]
            continue
    with open('./clustering/multigene_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return orf_increment, single_gene_cluster_complete, novel_genes_closest_refs


def add_gene_names_to_gbk(generics, gbk_dir):
    """
    Updates the Genbank files with the new gene names based on clustering

    :param generics: dict of isolates and gene names to update
    :param gbk_dir: str directory of Genbanks to update
    :return: None
    """
    logger = logging.getLogger('UpdateGenbank')
    gbk_files = os.listdir(gbk_dir)
    isolates = [os.path.splitext(isolate_id)[0] for isolate_id in gbk_files if isolate_id.endswith('.gbk')]
    for isolate in generics.keys():
        if isolate not in isolates:
            logger.error('Isolate ' + isolate + ' absent in ' + gbk_dir)
        else:
            genbank_file = gbk_dir + '/' + isolate + '.gbk'
            isolate_records = list(SeqIO.parse(genbank_file, 'genbank'))
            update_orf_dict = generics[isolate]
            locus_to_update = list(update_orf_dict.keys())
            modified_locus = []
            for rec in isolate_records:
                for feature in rec.features:
                    if feature.type != 'CDS':
                        continue
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    gene_name = feature.qualifiers['gene'][0]
                    try:
                        gene_synonym = feature.qualifiers['gene_synonym']
                    except KeyError:
                        gene_synonym = []
                    if locus_tag in locus_to_update and \
                            locus_tag not in modified_locus:
                        if designator.is_raw_ltag(gene_name):
                            feature.qualifiers['gene'][0] = update_orf_dict[locus_tag]['name']
                            if update_orf_dict[locus_tag]['pseudo']:
                                feature.qualifiers['pseudo'] = ['']
                                feature.qualifiers.pop('translation', None)
                        elif gene_name == update_orf_dict[locus_tag]['name']:
                            modified_locus.append(locus_tag)
                            continue
                        else:
                            logger.debug('Discordant assignment of gene name')
                            logger.debug('Original gene name: ' + gene_name)
                            logger.debug('New gene name: ' + update_orf_dict[locus_tag]['name'])
                            designator.append_qualifier(
                                feature.qualifiers,
                                'gene_synonym',
                                update_orf_dict[locus_tag]['name']
                            )
                        modified_locus.append(locus_tag)
                    elif locus_tag in locus_to_update and locus_tag in modified_locus:
                        logger.debug('Updated locus tag previously')
                    else:
                        continue
            if len(set(locus_to_update).intersection(set(modified_locus))) < len(locus_to_update):
                logger.warning('The following locus_tags are missing in the genbank file')
                logger.warning(set(locus_to_update).difference(set(modified_locus)))
            SeqIO.write(isolate_records, genbank_file, 'genbank')


def parseClustersUpdateGBKs(target_gffs, clusters, genomes_to_annotate, seq_ident, seq_covg):
    """
    Executes all functions to parse the clustering file and update
    all Genbanks

    :param gffs: str directory of GFFs
    :param genome_seqs: dict of sample name to SeqRecord genome sequence
    :param clusters: str cluster file
    :return: None
    """
    logger = logging.getLogger('ParseClusters')
    hybran_tmp_dir = config.hybran_tmp_dir
    global isolate_update_dictionary, isolate_sequences
    isolate_update_dictionary = {}
    logger.debug('Retrieving unique protein sequences from GFFs')
    # Run CD-HIT on cdss_protein-all.fasta as part of calling ref_seqs
    unique_protein_fastas, isolate_sequences = unique_seqs(annotations=target_gffs)
    logger.debug('Identifying reference and reference-transferred proteins')
    reference_protein_fastas = os.path.join(
        hybran_tmp_dir,
        'clustering',
        'ref_cdss_protein.fasta'
    )
    extractor.subset_fasta(
        inseq = unique_protein_fastas,
        outseq = reference_protein_fastas,
        match = designator.is_reference,
    )
    logger.debug('Identifying unannotated proteins (signified by absence of a gene name)')
    # Identify all the candidate novel genes in cdss_protein.fasta
    generic_genes_fp = os.path.join(hybran_tmp_dir,
                                'clustering',
                                'unannotated_seqs.fasta')
    extractor.subset_fasta(
        inseq = unique_protein_fastas,
        outseq = generic_genes_fp,
        match = designator.is_unannotated,
    )
    orf_increment = designator.find_next_increment(fasta=generic_genes_fp)
    logger.info('Parsing ' + clusters)
    clusters = parse_clustered_proteins(clustered_proteins=clusters,
                                        annotations=target_gffs)
    multi_gene_cluster = clusters[0]
    single_gene_cluster = clusters[1]
    candidate_novel_gene_cluster = clusters[2]
    unique_gene_cluster = clusters[3]
    candidate_novel_gene_cluster_complete = candidate_novel_gene_cluster.copy()
    single_gene_cluster_complete = single_gene_cluster.copy()
    logger.debug('Annotating genes that do not have a gene name but cluster with genes that have many gene names')
    orf_increment, \
        updated_single_gene_clusters, \
        novelty_report_multiref = multigene_clusters(in_dict=multi_gene_cluster,
                                                          single_gene_cluster_complete=single_gene_cluster_complete,
                                                          unannotated_fasta=generic_genes_fp,
                                                          orf_increment=orf_increment,
                                                          seq_ident=seq_ident,
                                                          seq_covg=seq_covg)
    logger.debug('Annotating genes that do not have a gene name but cluster with genes that have the same gene name')
    single_gene_clusters(updated_single_gene_clusters)
    orf_increment, novelty_report_noref = only_ltag_clusters(in_dict=candidate_novel_gene_cluster_complete,
                                       unannotated_fasta=generic_genes_fp,
                                       orf_increment=orf_increment,
                                       reference_fasta=reference_protein_fastas,
                                       seq_ident=seq_ident,
                                       seq_covg=seq_covg)
    logger.debug('Annotating genes that do not have a gene name and exist in a single isolate')
    orf_increment, novelty_report_singleton = singleton_clusters(singleton_dict=unique_gene_cluster,
                                       reference_fasta=reference_protein_fastas,
                                       unannotated_fasta=generic_genes_fp,
                                       orf_increment=orf_increment,
                                       seq_ident=seq_ident,
                                       seq_covg=seq_covg)
    with open('clustering/novelty_report.tsv','w') as novelty_report:
        print('\t'.join(['cluster_type',
                         'candidate_novel_gene',
                         'nearest_ref_match',
                         'metric',
                         'pct_aa_ident',
                         'pct_sub_covg',
                         'pct_qry_covg',
                         ]),
              file=novelty_report)
        for report in novelty_report_singleton + novelty_report_noref + novelty_report_multiref:
            print(report, file=novelty_report)
    logger.info('Updating Genbank files in ' + os.getcwd())
    add_gene_names_to_gbk(generics=isolate_update_dictionary,
                          gbk_dir=os.getcwd())
    logger.info('Preparing FASTA for eggNOG functional assignments')
    with open(os.path.join(hybran_tmp_dir,'eggnog_seqs.fasta'), 'w') as eggnog_seqs:
        prepare_for_eggnog(unannotated_seqs=generic_genes_fp,
                           outfile = eggnog_seqs)
