
import os
import re
import tempfile
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastpCommandline
from . import config


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
    underscore_re = re.compile('_[0-9]$')

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
                            if len(column) >= 8 and 'CDS' in column[2]:
                                info = column[8].split(';')
                                gene_id = ''.join([i.split('=')[1] for i in info if i.startswith('ID=')])
                                locus_tag = ''.join([i.split('=')[1] for i in info if i.startswith('locus_tag=')])
                                gene = ','.join([i.split('=')[1] for i in info if i.startswith('gene=')])
                                if not locus_tag.startswith('Rv') and not locus_tag.startswith('L'):
                                    gene = locus_tag
                                isolate_id_ltag[gene_id] = (locus_tag, gene)
                        gff_dictionary[isolate_id] = isolate_id_ltag
            except IOError:
                continue
        return gff_dictionary

    # Different gene names/Rvs in cluster
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
            isolates_ids = line.rstrip('\n').split('\t')
            representative_seq_id = isolates_ids[0].split(': ')[1]
            rep_seq_id = representative_seq_id.split('-')[-1]
            rep_isolate = representative_seq_id.replace('-' + rep_seq_id, '')
            representative_ltag_gene_tup = gffs[rep_isolate][rep_seq_id]
            representative = ','.join([rep_isolate, ','.join(representative_ltag_gene_tup)])

            representative_fasta_list.append([rep_isolate] + list(representative_ltag_gene_tup))
            gene_cluster[representative] = []

            cluster_list = []
            cluster_list_w_isolate = []
            cluster_list.append(representative_ltag_gene_tup)
            cluster_list_w_isolate.append([rep_isolate] + list(representative_ltag_gene_tup))
            for isolate_gene_id in isolates_ids[1:]:
                gene_id = isolate_gene_id.split('-')[-1]
                isolate = isolate_gene_id.replace('-' + gene_id, '')
                cluster_list.append(gffs[isolate][gene_id])
                cluster_list_w_isolate.append([isolate] + list(gffs[isolate][gene_id]))
                gene_cluster[representative].append(','.join([isolate, ','.join(gffs[isolate][gene_id])]))
            number_genes, number_locus_tags = different_genes(cluster_list)
            unique_genes = []
            unique_locus = []
            for data in number_genes:
                if not data.startswith('L'):
                    unique_genes.append(data)
            for data in number_locus_tags:
                if not data.startswith('L'):
                    unique_locus.append(data)
            # Getting clusters with only 1 gene that is an L tag or an underscore
            if (len(unique_genes) == 1 or len(unique_locus) == 1) and \
                    any([gene[1].startswith('L') for gene in cluster_list]):
                same_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            elif (len(isolates_ids) == 1 and
                    (representative_ltag_gene_tup[1].startswith('L')) or
                     underscore_re.search(representative_ltag_gene_tup[1])):
                unique_genes_list.append([rep_isolate] + list(representative_ltag_gene_tup))
            # All different with L tags
            elif (len(unique_genes) > 1 or len(unique_locus) > 1) and \
                    any([gene[1].startswith('L') for gene in cluster_list]):
                different_genes_cluster_w_ltags[representative] = cluster_list_w_isolate
            #
            # L tag only clusters
            elif all(locus[1].startswith('L') for locus in cluster_list):
                l_tag_only_clusters[representative] = cluster_list_w_isolate
    return [different_genes_cluster_w_ltags, same_genes_cluster_w_ltags, l_tag_only_clusters,
            unique_genes_list]


def get_top_hit(all_hits_dict):
    """
    Gets best BLASTP hit based on amino acid sequence identity

    :param all_hits_dict: This dictionary consists of all Blast hits
    with corresponding unannotated locus_tag that have
    passed the defined thresholds for amino acid identity and coverage
    :return: list of only the elements of top hits with the locus tag.
    The key is the 'L(2)_' locus tag and the value is
     corresponding top hit in H37Rv. If multiple top hits are
     present with same identity and coverage values, the value
     is all the top Rv hits separated by ':'
    """
    top_hit_identity = 0
    for gene in all_hits_dict.keys():
        if all_hits_dict[gene] > top_hit_identity:
            top_hit_identity = all_hits_dict[gene]
            top_hit = gene
    return top_hit


def identify_top_hits(blast_output_file, identity, coverage):
    """
    Iterates over BLASTP output and determines all the best
    hits based on identity and alignment coverage

    :param blast_output_file: Output file from Blastp
    :return: Dictionary of top hits for each unannotated locus that
    passes the identity and coverage threshold
             If none exist, a list of three report strings for the next-closest matches by each category are reported
    """
    all_hits_dict = {}
    subcritical_hits_iden = {}
    subcritical_hits_scov = {}
    subcritical_hits_qcov = {}
    blast_output_raw = blast_output_file.split('\n')
    rv_hit = False
    for line in blast_output_raw:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        line_elements = line.split('\t')
        if line_elements and line_elements != [' ']:
            iden = float(line_elements[5])
            qcov = (float(line_elements[4]) / float(line_elements[1])) * 100.0
            scov = (float(line_elements[4]) / float(line_elements[3])) * 100.0
            if not line_elements[2].startswith('L'):
                if '|' in line_elements[2]:
                    [seqid, locus_tag, gene_name] = line_elements[2].split('|')
                    if gene_name:
                        corresponding_rv_hit = gene_name
                    else:
                        corresponding_rv_hit = locus_tag
                else:
                    corresponding_rv_hit = line_elements[2]

                if iden >= identity and qcov >= coverage and scov >= coverage:
                    rv_hit = True
                    all_hits_dict[corresponding_rv_hit] = iden
                # Rv hits that don't meet the threshold.
                # Keep track of them so we know how close we were to getting a good match
                else:
                    subcritical_hits_iden[corresponding_rv_hit] = iden
                    subcritical_hits_scov[corresponding_rv_hit] = scov
                    subcritical_hits_qcov[corresponding_rv_hit] = qcov
            else:
                continue
    if rv_hit:
        top_hit = get_top_hit(all_hits_dict)
        best_subcriticals = []
    else:
        top_hit = None
        if subcritical_hits_iden:
            best_subc_iden = max(subcritical_hits_iden, key=lambda _: subcritical_hits_iden[_])
            best_subc_scov = max(subcritical_hits_scov, key=lambda _: subcritical_hits_scov[_])
            best_subc_qcov = max(subcritical_hits_qcov, key=lambda _: subcritical_hits_qcov[_])
            best_subcriticals = [
                '\t'.join([_[0],
                           _[1],
                           str(subcritical_hits_iden[_[0]]),
                           str(subcritical_hits_scov[_[0]]),
                           str(subcritical_hits_qcov[_[0]])])
                for _ in zip([best_subc_iden, best_subc_scov, best_subc_qcov],
                             ['identity','scov','qcov'])]
        else:
            best_subcriticals = []
    return top_hit, best_subcriticals, all_hits_dict


def update_dictionary_ltag_assignments(isolate_id, isolate_ltag, new_gene_name):
    """
    Updates a global dictionary with new gene name assignments

    :param isolate_id: str isolate ID
    :param isolate_ltag: str locus tag to update
    :param new_gene_name: str gene name to assign
    :return: None
    """
    logger = logging.getLogger('NameGene')
    if isolate_id not in isolate_update_dictionary.keys():
        isolate_update_dictionary[isolate_id] = {}
        isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
        logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
    else:
        isolate_dict_added = isolate_update_dictionary[isolate_id]
        if isolate_ltag not in isolate_dict_added.keys():
            logger.debug(isolate_ltag + ' in ' + isolate_id + ' becomes ' + new_gene_name)
            isolate_update_dictionary[isolate_id][isolate_ltag] = new_gene_name
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
    if not (rep_locus.startswith('L') and rep_gene_name.startswith('L')):
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
        gene_isolate_id = gene[0]
        gene_locus = gene[1]
        gene_name = gene[2]
        if gene_locus.startswith('L') and gene_name.startswith('L'):
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
        if not (gene_locus.startswith('L') and gene_name.startswith('L')):
            annotated += 1
    if annotated == len(cluster_list):
        return True
    else:
        return False


def find_largest_mtb_increment(unannotated_fasta):
    """
    Based on a given unannotated FASTA, identifies the highest
    numbered increment of MTBs

    :param unannotated_fasta: FASTA file name
    :return: int
    """
    mtbs = []
    for record in SeqIO.parse(unannotated_fasta, 'fasta'):
        mtbs.append(record.id)
    if mtbs:
        largest_mtb = sorted(mtbs)[-1]
        return int(largest_mtb.replace('MTB', ''))
    else:
        return 1


def ref_seqs(gbk_dir):
    """
    Creates a FASTA with unique sequences based on a given directory
    of reference Genbank files. Uses GFFs from referene and input
    genomes to create content for the output ref_cdss_protein-all.fasta

    :param gbk_dir: str directory name
    :return: str FASTA file name and a dictionary of dictionaries
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    def create_allseq_dict(fa):
        """
        Parses FASTA input, creates dictionary with
        sequence ID as key and the sequence as the value

        :param fa: FASTA file name
        :return: dict of record IDs and sequences
        """
        fasta_dict = {}
        fasta_descrip = {}
        for record in SeqIO.parse(fa, "fasta"):
            sequence = str(record.seq)
            fasta_dict[record.id] = sequence
            fasta_descrip[record.id] = record.description
        return fasta_dict, fasta_descrip

    def run_cdhit(nproc, input, output):
        """
        Runs cdhit with a sequence identity threshold of 0.99.

        CDHIT outputs a fasta file of representative sequences (ref_cdss_protein.fasta)
        and a text file of list of clusters (ref_cdss_protein.fasta.clstr).

        :param nproc: str number of processors to use
        :param input: str input file name
        :param output: str output file name
        :return: subprocess.stdout
        """
        cmd = ['cd-hit',
               '-i', input,
               '-o', output,
               '-c', '0.99',
               '-T', str(nproc),
               '-g', '1',
               '-s', '1.0',
               '-d', '256',
               '-A', '1.0']
        out = subprocess.run(cmd, stdout=subprocess.PIPE)
        return out

    def create_reps_dict(in_clusters):
        """
        Creates a dictionary with the representative
        sequence ID from each cluster as the key
        and the corresponding sequences as the value.

        :param  in_clusters: str file name of clusters
        :return: dict of representatives and the CDHIT clusters and a list of representatives
        """
        clusters = []
        rep_dict = {}
        rep = ''
        reps = []
        with open(in_clusters, 'r') as f:
            for line in f:
                # At new cluster
                if line.startswith('>'):
                    if not rep:
                        continue
                    rep_dict[rep] = clusters
                    clusters = []
                    continue
                # Representatives of cluster
                else:
                    if line.rstrip('\n').endswith('*'):
                        rep = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                        reps.append(rep)
                    else:
                        gene = line.rstrip('\n').split()[2].replace('>', '').replace('...', '')
                        clusters.append(gene)
        rep_dict[rep] = clusters
        return rep_dict, reps

    def create_reps_fasta(output, reps, rep_seq_id_dict, description_dict):
        """
        Creates a FASTA file that has sequence ID
        and sequence for each representative

        :param output: str output file name
        :param reps: list of representatives
        :param rep_seq_id_dict: dict of representatives and sequences
        :return: None
        """
        seq_list = []
        for id_key in reps:
            record = SeqRecord(Seq(rep_seq_id_dict[id_key]), id=id_key, description=description_dict[id_key])
            seq_list.append(record)
        SeqIO.write(seq_list, output, "fasta")

    def cd_hit(nproc, fasta, out):
        """
        Executes CDHIT

        :param nproc: str number of processors
        :param fasta: str FASTA file name
        :param out: str output file name
        :return: None
        """
        logger = logging.getLogger('CDHIT')
        logger.debug('Running CDHIT on reference annotations')
        cdhit_stdout = run_cdhit(nproc, fasta, out)
        OGdict, description_dict = create_allseq_dict(fasta)
        REPdict, rep_list = create_reps_dict(out + ".clstr")
        create_reps_fasta(out, rep_list, OGdict, description_dict)

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
    for gff in gbk_dir:
        if gff.endswith('.gff'):
            raw_out = grep_seqs(gff)
            gff_name = os.path.splitext(os.path.basename(gff))[0]
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
    all_proteins = hybran_tmp_dir + '/clustering/ref_cdss_protein-all.fasta'
    final_proteins = hybran_tmp_dir + '/clustering/ref_cdss_protein.fasta'
    # Write ref_cdss_protein-all.fasta
    with open(all_proteins, 'w') as ref_cds_out:
        for s in protein_cds:
            SeqIO.write(s, ref_cds_out, 'fasta')
    cd_hit(1, all_proteins, final_proteins)
    return final_proteins, isolate_seqs


def blast(subject, stdin_seq):
    """
    BLASTP the stdin_seq against the subject file

    :param subject: str FASTA file name
    :param stdin_seq: str amino acid sequence
    :return: list BLAST output
    """
    blast_command = NcbiblastpCommandline(subject=subject,
                                          outfmt='"7 qseqid qlen sseqid slen qlen length pident qcovs"')
    stdout, stderr = blast_command(stdin=str(stdin_seq))
    return stdout


def find_unannotated_genes(reference_protein_fasta):
    """
    Identifies all unannotated sequences in the reference proteome by
    looking for all the MTB#### genes in the given FASTA.

    :param reference_protein_fasta: str FASTA file
    :return: str unannotated FASTA file name
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    unannotated_seqs = []
    for record in SeqIO.parse(reference_protein_fasta, 'fasta'):
        if record.id.startswith('MTB'):
            unannotated_seqs.append(record)
    with open(hybran_tmp_dir + '/clustering/unannotated_seqs.fasta', 'w') as \
            out:
        for s in unannotated_seqs:
            SeqIO.write(s, out, 'fasta')
    return hybran_tmp_dir + '/clustering/unannotated_seqs.fasta'


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


def singleton_clusters(singleton_dict, reference_fasta, unannotated_fasta, mtb_increment, seq_ident, seq_covg):
    """
    If a cluster has only one gene (with no gene names),
    then BLAST representative sequence to H37Rv and if there
    is a hit with specified amino acid and coverage thresholds,
    all candidate novel genes in the cluster is annotated
    with the H37Rv gene. If the representative does not hit a H37Rv,
    assign a MTB locus tag to the genes in the  cluster.

    :param singleton_dict: dict of clusters
    :param reference_fasta: str reference FASTA proteome
    :param unannotated_fasta: str unannotated FASTA file
    :param mtb_increment: int increment to start new MTBs
    :return: int updated increment of MTBs
    """
    logger = logging.getLogger('SingletonClusters')

    logger.debug('Number of singleton clusters with single genes: ' + str(len(singleton_dict)))

    new_unannotated_seqs = []
    out_list = []
    novel_genes_closest_refs = []
    for single_gene in singleton_dict:
        isolate_id = single_gene[0]
        locus_tag = single_gene[1]
        gene_name = single_gene[2]
        out_list.append(str(single_gene))
        if gene_name.startswith('L') and locus_tag.startswith('L'):
            gene_sequence = isolate_sequences[isolate_id][locus_tag]
            stdout = blast(reference_fasta, gene_sequence)
            top_hit, best_subcriticals, all_hits = identify_top_hits(stdout, identity=seq_ident, coverage=seq_covg)
            assign_mtb = False
            name_to_assign = top_hit
            if top_hit is None:
                assign_mtb = True
                if os.stat(unannotated_fasta).st_size:
                    stdout_2 = blast(unannotated_fasta, gene_sequence)
                    top_hit_mtb, best_subcriticals_mtb, all_hits_mtb = identify_top_hits(stdout_2, identity=seq_ident, coverage=seq_covg)
                    best_subcriticals += best_subcriticals_mtb
                    if top_hit_mtb is not None:
                        assign_mtb = False
                        name_to_assign = top_hit_mtb
            if assign_mtb:
                mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                mtb_increment = mtb_increment + 1
                seq_record = SeqRecord(gene_sequence, id=mtb_id, description='False')
                new_unannotated_seqs.append(seq_record)
                name_to_assign = mtb_id
                best_subcriticals = ['\t'.join(['singleton', mtb_id, hit])
                                     for hit in best_subcriticals]
                novel_genes_closest_refs += best_subcriticals
                if not novel_genes_closest_refs:
                    novel_genes_closest_refs = ['\t'.join(['singleton', mtb_id])]
            update_dictionary_ltag_assignments(isolate_id, locus_tag, name_to_assign)
        else:
            continue
    with open(unannotated_fasta, 'a') as mtb_fasta:
        for s in new_unannotated_seqs:
            SeqIO.write(s, mtb_fasta, 'fasta')
    with open('./clustering/singleton_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return mtb_increment, novel_genes_closest_refs


def only_ltag_clusters(in_dict, reference_fasta, unannotated_fasta, mtb_increment, seq_ident, seq_covg):
    """
    If a cluster has only candidate novel genes (with no gene names),
    then BLAST representative sequence to H37Rv
    and if there is a hit with specified amino acid
    and coverage thresholds, all candidate novel genes in the cluster
    is annotated with the H37Rv gene. If the representative does
    not hit a H37Rv, assign a MTB locus tag to the genes
    in the cluster.

    :param in_dict: dict of clusters
    :param reference_fasta: str reference FASTA proteome
    :param unannotated_fasta: str unannotated FASTA file
    :param mtb_increment: int increment to start new MTBs
    :return: int updated increment of MTBs
    """
    hybran_tmp_dir = config.hybran_tmp_dir
    logger = logging.getLogger('NewGeneClusters')

    logger.debug('Number of clusters with no gene names: ' + str(len(list(in_dict.keys()))))

    new_unannotated_genes = []
    rep_records = []
    out_list = []
    novel_genes_closest_refs = []
    for rep_gene in in_dict.keys():
        rep_isolate_id = rep_gene.split(',')[0]
        rep_locus = rep_gene.split(',')[1]
        rep_gene_name = rep_gene.split(',')[2]
        rep_sequence = isolate_sequences[rep_isolate_id][rep_locus]
        rep_record = SeqRecord(rep_sequence, id='|'.join([rep_isolate_id, rep_locus, rep_gene_name]),
                               description='')
        rep_records.append(rep_record)
        out_list.append(rep_gene + ': ' + str(in_dict[rep_gene]))

    rep_fp = tempfile.NamedTemporaryFile(suffix='.fasta',
                                         dir=hybran_tmp_dir,
                                         delete=False, mode='w')
    rep_temp_fasta = rep_fp.name
    with open(rep_temp_fasta, 'w') as fasta_tmp:
        for s in rep_records:
            SeqIO.write(s, fasta_tmp, 'fasta')
    blastcmd = NcbiblastpCommandline(subject=reference_fasta,
                                     query=rep_temp_fasta,
                                     outfmt='"6 qseqid qlen sseqid slen qlen length pident qcovs"')
    stdout, stderr = blastcmd()
    rep_fp.close()
    hits = {}
    for line in stdout.split('\n'):
        if line:
            try:
                hits[line.split('\t')[0]] += [line]
            except KeyError:
                hits[line.split('\t')[0]] = [line]
    for qid, lines in hits.items():
        top_hit, best_subcriticals, all_hits = identify_top_hits('\n'.join(lines), identity=seq_ident, coverage=seq_covg)
        isolate = qid.split('|')[0]
        locus = qid.split('|')[1]
        gene = qid.split('|')[2]
        assign_mtb = False
        name_to_assign = top_hit
        if top_hit is None:
            assign_mtb = True
            if os.stat(unannotated_fasta).st_size:
                rep_seq = isolate_sequences[qid.split('\t')[0].split('|')[0]][qid.split('\t')[0].split('|')[1]]
                stdout_2 = blast(unannotated_fasta, rep_seq)
                top_hit_mtb, best_subcriticals_mtb, all_hits_mtb = identify_top_hits(stdout_2, identity=seq_ident, coverage=seq_covg)
                best_subcriticals += best_subcriticals_mtb
                if top_hit_mtb is not None:
                    assign_mtb = False
                    name_to_assign = top_hit_mtb
        if assign_mtb:
            mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
            mtb_increment = mtb_increment + 1
            rep_sequence = isolate_sequences[qid.split('\t')[0].split('|')[0]][qid.split('\t')[0].split('|')[1]]
            seq_record = SeqRecord(rep_sequence, id=mtb_id, description='False')
            new_unannotated_genes.append(seq_record)
            name_to_assign = mtb_id
            best_subcriticals = ['\t'.join(['no_refs', mtb_id, hit])
                                 for hit in best_subcriticals]
            novel_genes_closest_refs += best_subcriticals
            if not novel_genes_closest_refs:
                novel_genes_closest_refs = ['\t'.join(['no_refs', mtb_id])]
        update_dictionary_ltag_assignments(isolate, locus, name_to_assign)
        for other_genes_in_cluster in in_dict[','.join([isolate, locus, gene])]:
            update_dictionary_ltag_assignments(other_genes_in_cluster[0], other_genes_in_cluster[1], name_to_assign)
    with open(unannotated_fasta, 'a') as mtb_fasta:
        for s in new_unannotated_genes:
            SeqIO.write(s, mtb_fasta, 'fasta')
    with open('./clustering/onlyltag_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return mtb_increment, novel_genes_closest_refs


def single_gene_clusters(single_gene_dict):
    """
    Handling clusters with only one Rv/gene name in cluster
    If cluster has candidate novel genes (L_***** or L2_*****)
    clustered together with a H37Rv annotation, update
    all candidate novel genes in this cluster with the H37Rv gene name.

    :param single_gene_dict: dict of clusters
    :return: None
    """
    logger = logging.getLogger('ClustersWithOneGeneName')

    rep_ltag_keys = []
    logger.debug('Number of clusters with single genes: ' + str(len(list(single_gene_dict.keys()))))

    for key_sgc in single_gene_dict:
        rep_ltag = False
        key_elements_sgc = key_sgc.split(',')
        if key_elements_sgc[1].startswith('Rv'):
            if key_elements_sgc[2]:
                gene_to_add = key_elements_sgc[2]
            else:
                gene_to_add = key_elements_sgc[1]
        elif key_elements_sgc[1].startswith('L') and key_elements_sgc[2].startswith('L'):
            rep_ltag = True
            rep_ltag_keys.append(key_sgc)
        else:
            gene_to_add = key_elements_sgc[2]
        # If representative sequence is an Rv and not an L-tag
        if not rep_ltag:
            genes_in_cluster = single_gene_dict[key_sgc]
            for gene in genes_in_cluster:
                if gene[1].startswith('L') and gene[2].startswith('L'):
                    update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                else:
                    continue
        # If representative sequence is a L-tag
        else:
            genes_in_cluster = single_gene_dict[key_sgc]
            gene_to_add = ''
            for gene in genes_in_cluster:
                if gene[1].startswith('L') and gene[2].startswith('L'):
                    continue
                else:
                    if not gene[1].startswith('L'):
                        if gene[2]:
                            gene_to_add = gene[2]
                        else:
                            gene_to_add = gene[1]
                    break
            if len(gene_to_add) != 0:
                for gene in genes_in_cluster:
                    if gene[1].startswith('L') and gene[2].startswith('L'):
                        update_dictionary_ltag_assignments(gene[0], gene[1], gene_to_add)
                    else:
                        continue


def multigene_clusters(in_dict, single_gene_cluster_complete, unannotated_fasta, mtb_increment, seq_ident, seq_covg):
    """
    If a cluster has multiple H37Rv genes and L_tags,
    BLAST the L_tags to all genes in the cluster that are H37Rv
    and annotate L_tag with the top hit.

    :param in_dict: dict of clusters
    :param single_gene_cluster_complete: dict of clusters that have a single sequence
    :param unannotated_fasta: str unannotated FASTA file
    :param mtb_increment: int increment to start new MTBs
    :return: int updated increment of MTBs
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
            if gene_in_cluster[1] in gene_elements or gene_in_cluster[2] in gene_elements:
                continue
            elif gene_in_cluster[1].startswith('L') and gene_in_cluster[2].startswith('L'):
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
                    unannotated_gene_seq = isolate_sequences[unannotated_gene_isolate][unannotated_gene_locus]
                    stdout = blast(fasta_fp_to_blast, unannotated_gene_seq)
                    top_hit, best_subcriticals, all_hits = identify_top_hits(stdout, identity=seq_ident, coverage=seq_covg)
                    assign_mtb = False
                    name_to_assign = top_hit
                    if top_hit is None:
                        assign_mtb = True
                        if os.stat(unannotated_fasta).st_size:
                            stdout_2 = blast(unannotated_fasta, unannotated_gene_seq)
                            top_hit_mtb, best_subcriticals_mtb, all_hits_mtb = identify_top_hits(stdout_2, identity=seq_ident, coverage=seq_covg)
                            best_subcriticals += best_subcriticals_mtb
                            if top_hit_mtb is not None:
                                assign_mtb = False
                                name_to_assign = top_hit_mtb
                    if assign_mtb:
                        mtb_id = 'MTB' + "%04g" % (int('0001') + mtb_increment)
                        mtb_increment = mtb_increment + 1
                        seq_record = SeqRecord(unannotated_gene_seq, id=mtb_id, description='False')
                        with open(unannotated_fasta, 'a') as mtb_fasta:
                            SeqIO.write(seq_record, mtb_fasta, 'fasta')
                        name_to_assign = mtb_id
                        best_subcriticals = ['\t'.join(['multiref', mtb_id, hit])
                                             for hit in best_subcriticals]
                        novel_genes_closest_refs += best_subcriticals
                        if not novel_genes_closest_refs:
                            novel_genes_closest_refs = ['\t'.join(['multiref', mtb_id])]
                    update_dictionary_ltag_assignments(unannotated_gene_isolate, unannotated_gene_locus, name_to_assign)
        else:
            if len(unassigned_l_tags) > 0:
                single_gene_cluster_complete[gene] = in_dict[gene]
            continue
    with open('./clustering/multigene_clusters.txt', 'w') as f:
        for line in out_list:
            f.write(str(line) + '\n')
    return mtb_increment, single_gene_cluster_complete, novel_genes_closest_refs


def add_gene_names_to_gbk(mtbs, gbk_dir):
    """
    Updates the Genbank files with the new gene names based on clustering

    :param mtbs: dict of isolates and gene names to update
    :param gbk_dir: str directory of Genbanks to update
    :return: None
    """
    logger = logging.getLogger('UpdateGenbank')
    gbk_files = os.listdir(gbk_dir)
    isolates = [os.path.splitext(isolate_id)[0] for isolate_id in gbk_files if isolate_id.endswith('.gbk')]
    for isolate in mtbs.keys():
        if isolate not in isolates:
            logger.error('Isolate ' + isolate + ' absent in ' + gbk_dir)
        else:
            genbank_file = gbk_dir + '/' + isolate + '.gbk'
            isolate_records = list(SeqIO.parse(genbank_file, 'genbank'))
            update_mtb_dict = mtbs[isolate]
            locus_to_update = list(update_mtb_dict.keys())
            modified_locus = []
            rec_num = 1
            for rec in isolate_records:
                if rec_num > 1:
                    break
                rec_num += 1
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
                        if gene_name.startswith('L'):
                            feature.qualifiers['gene'] = update_mtb_dict[locus_tag]
                        elif gene_name == update_mtb_dict[locus_tag]:
                            modified_locus.append(locus_tag)
                            continue
                        else:
                            logger.debug('Discordant assignment of gene name')
                            logger.debug('Original gene name: ' + gene_name)
                            logger.debug('New gene name: ' + update_mtb_dict[locus_tag])
                            if 'gene_synonym' in feature.qualifiers.keys():
                                gene_synonym.append(update_mtb_dict[locus_tag])
                            else:
                                feature.qualifiers['gene_synonym'] = [update_mtb_dict[locus_tag]]
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
    :param clusters: str cluster file
    :return: None
    """
    logger = logging.getLogger('ParseClusters')
    hybran_tmp_dir = config.hybran_tmp_dir
    global isolate_update_dictionary, isolate_sequences
    isolate_update_dictionary = {}
    logger.debug('Retrieving reference protein sequences from GFFs')
    # Run CD-HIT on ref_cdss_protein-all.fasta as part of calling ref_seqs
    reference_protein_fastas, isolate_sequences = ref_seqs(gbk_dir=target_gffs)
    logger.debug('Identifying unannotated proteins (signified by absence of a gene name)')
    # Identify all the MTB#### genes in ref_cdss_protein-all.fasta
    mtb_genes_fp = find_unannotated_genes(reference_protein_fasta=reference_protein_fastas)
    mtb_increment = find_largest_mtb_increment(unannotated_fasta=mtb_genes_fp)
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
    mtb_increment, \
        updated_single_gene_clusters, \
        novelty_report_multiref = multigene_clusters(in_dict=multi_gene_cluster,
                                                          single_gene_cluster_complete=single_gene_cluster_complete,
                                                          unannotated_fasta=mtb_genes_fp,
                                                          mtb_increment=mtb_increment,
                                                          seq_ident=seq_ident,
                                                          seq_covg=seq_covg)
    logger.debug('Annotating genes that do not have a gene name but cluster with genes that have the same gene name')
    single_gene_clusters(updated_single_gene_clusters)
    mtb_increment, novelty_report_noref = only_ltag_clusters(in_dict=candidate_novel_gene_cluster_complete,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment,
                                       reference_fasta=reference_protein_fastas,
                                       seq_ident=seq_ident,
                                       seq_covg=seq_covg)
    logger.debug('Annotating genes that do not have a gene name and exist in a single isolate')
    mtb_increment, novelty_report_singleton = singleton_clusters(singleton_dict=unique_gene_cluster,
                                       reference_fasta=reference_protein_fastas,
                                       unannotated_fasta=mtb_genes_fp,
                                       mtb_increment=mtb_increment,
                                       seq_ident=seq_ident,
                                       seq_covg=seq_covg)
    with open('clustering/novelty_report.tsv','w') as novelty_report:
        print('\t'.join(['cluster_type',
                         'candidate_novel_gene',
                         'nearest_ref_match',
                         'metric',
                         '% AA sequence identity',
                         '% subject(ref) coverage',
                         '% query coverage',
                         ]),
              file=novelty_report)
        for report in novelty_report_singleton + novelty_report_noref + novelty_report_multiref:
            print(report, file=novelty_report)
    logger.info('Updating Genbank files in ' + os.getcwd())
    add_gene_names_to_gbk(mtbs=isolate_update_dictionary,
                          gbk_dir=os.getcwd())
    logger.info('Preparing FASTA for eggNOG functional assignments')
    with open(os.path.join(hybran_tmp_dir,'eggnog_seqs.fasta'), 'w') as eggnog_seqs:
        prepare_for_eggnog(unannotated_seqs=mtb_genes_fp,
                           outfile = eggnog_seqs)