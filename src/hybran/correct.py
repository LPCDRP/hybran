from collections import defaultdict
from copy import deepcopy
from glob import glob
import itertools
from multiprocessing import Pool
import os
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord
from frozendict import frozendict
import networkx as nx

from . import (
    annomerge,
    converter,
    BLAST,
    designator,
    extractor,
    fileManager,
)
from .annomerge import merge, thunderdome
from .bio import (
    AutarkicSeqFeature,
    sort_features,
    translate,
)
from .config import cnf
from .compare import compare, cross_examine
from .parseClustering import resolve_clusters
from .pseudoscan import pseudoscan
from .util import keydefaultdict
# from .util import mpbreakpoint


genes = defaultdict(dict)
gene_names = defaultdict(dict)
gene_seqs = defaultdict(dict)
strain_seqs = defaultdict(dict)
iso = []
logs_dir = None
seq_files = {}

dummy_bed_score = '0'
sign = lambda x: '+' if x==1 else '-'


def prepare_node(strain_id, identifier, feature):
    """
    Create a tuple for use in networkx Graph.add_node in order to consistently add features as nodes in the graph.
    """
    return (
        identifier,
        {
            'strain':strain_id,
            'annotation':feature,
        }
    )

def main(args):
    global genes
    global gene_names
    global gene_seqs
    global strain_seqs
    global iso
    global logs_dir
    global seq_files
    global ref_annotation

    cnf.blast.min_coverage = args.blast_min_coverage
    cnf.blast.min_identity = args.blast_min_identity

    try:
        os.makedirs(args.outdir, exist_ok=True)
    except:
        sys.exit(f"Could not create directory {args.outdir}")

    logs_dir = os.path.join(args.outdir, "logs")
    reports_dir = os.path.join(args.outdir, "reports")

    mcps_fn = os.path.join(reports_dir, "mcps.tsv")
    renames_fn = os.path.join(reports_dir, "renames.tsv")
    additions_fn = os.path.join(reports_dir, "additions.bed")

    n_inputs = len(args.annotations)
    input_type = None

    # Check the conditionally required arguments.
    if (
            n_inputs == 1
            and fileManager.is_hybran_output_dir(args.annotations[0])
    ):
        input_type = 'hybran_result'

        ann_files = glob.glob(
            os.path.join(args.annotations[0], '*.gbk')
        )

        if not args.references:
            args.references = glob.glob(os.path.join(
                args.annotations[0],
                'unified-refs',
                '*.gbk'
            ))
        else:
            args.references = fileManager.file_list(
                args.references,
                file_type='genbank',
            )
    else:
        if n_inputs == 1 and args.annotations[0].endswith('.bed'):
            input_type = 'bed'
            coords_file = args.annotations[0]
        else:
            input_type = 'genbank'
            ann_files = fileManager.file_list(args.annotations, file_type='genbank')

        missing_args = []
        if input_type=='bed' and not args.seq_dir:
            missing_args.append('-s/--seq-dir')
        if not args.references:
            missing_args.append('-r/--references')
        if missing_args:
            sys.exit(f"ERROR: missing arguments {', '.join(missing_args)}, which are required when not reading from a hybran output folder")

        args.references = fileManager.file_list(args.references, file_type='genbank')


    ref_annotation = extractor.load_gbks(
        args.references,
        feature_types=['CDS'],
    )
    annomerge.ref_annotation = ref_annotation

    strain_contig_records = {}

    if args.seq_dir:
        for fasta in fileManager.file_list([args.seq_dir], file_type='fasta'):
            strain_id = os.path.splitext(os.path.basename(fasta))[0]
            seq_files[strain_id] = fasta
            strain_contig_records[strain_id] = {
                f"{strain_id}.{record.id}":record for record in SeqIO.parse(fasta, "fasta")
            }
    # The only way we should get here is if our input is a Hybran result folder or genbank input.
    # We are assured to have the full sequences included in the genbank files.
    else:
        seq_dir = os.path.join(cnf.tmpdir, 'genomes')
        os.makedirs(seq_dir)
        for annotation in ann_files:
            curr_strain = os.path.splitext(os.path.basename(annotation))[0]
            strain_contig_records[curr_strain] = {
                record.id:record for record in SeqIO.parse(annotation, 'genbank')
            }
            curr_strain_fasta = os.path.join(seq_dir, f'{curr_strain}.fasta')
            _, strain_gene_seqs, _ = extractor.fastaFromGbk(
                annotation,
                out_cds=None,
                out_genome=curr_strain_fasta,
                identify=lambda f: f"{extractor.get_ltag(f)}|{extractor.get_gene(f)}|{f.location.strand}",
            )
            seq_files[curr_strain] = curr_strain_fasta
            gene_seqs[curr_strain] = {f.id.split('|')[0]:f for f in strain_gene_seqs}


    cnf.genetic_code = extractor.get_genetic_code(args.references[0])

    node_data = []
    unique_gene_names = set()

    if input_type == 'bed':
        sort_proc = subprocess.run([
            "sort",
            "-k1,1",
            "-k2,2n",
            coords_file,
        ], capture_output=True, check=True, text=True)
        last_strain = None
        cds_counter = defaultdict(lambda : 1)
        for line in sort_proc.stdout.strip().split('\n'):
            (
                strain_chrom,
                start,
                end,
                gene,
                _,
                strand_sign,
            ) = line.split('\t')
            strand = 1 if strand_sign == '+' else -1
            isolate_id, contig = os.path.splitext(strain_chrom)
            strain_contig_record = strain_contig_records[isolate_id][strain_chrom]
            # stand-in for locus tag
            cds_id = f"{isolate_id}_{cds_counter[isolate_id]}"
            cds_counter[isolate_id] += 1
            feature = AutarkicSeqFeature(
                type='CDS',
                location=SimpleLocation(
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    ref=strain_chrom,
                ),
                qualifiers={
                    'locus_tag':[cds_id],
                    'gene':gene,
                },
            )
            feature.references = {strain_chrom:strain_contig_record.seq}
            genes[isolate_id][cds_id] = feature
            gene_seqs[isolate_id][cds_id] = SeqRecord(
                translate(feature.extract(strain_contig_records[isolate_id][strain_chrom].seq), table=cnf.genetic_code),
                id=f"{cds_id}|{gene}|{strand}",
                description='',
            )
            gene_names[isolate_id][cds_id] = gene
            strain_contig_record.features.append(genes[isolate_id][cds_id])
            last_strain = isolate_id
            unique_gene_names.add(gene)

            node_data.append(prepare_node(isolate_id, cds_id, feature))
    else:
        for annotation in ann_files:
            isolate_id = os.path.splitext(os.path.basename(annotation))[0]
            for record in SeqIO.parse(annotation, "genbank"):
                for i, f in enumerate(record.features):
                    if f.type != 'CDS':
                        continue
                    f = AutarkicSeqFeature.fromSeqFeature(f)
                    for part in f.location.parts:
                        part.ref = record.id
                    f.references = {record.id: record.seq}
                    ltag = f.qualifiers['locus_tag'][0]
                    if 'inference' in f.qualifiers:
                        for inf_note in f.qualifiers['inference']:
                            if (
                                    inf_note.startswith('similar to ')
                                    and (
                                        inf_note.endswith(':RATT')
                                        or inf_note.endswith(':blastp')
                                    )
                            ):
                                ref_id = inf_note.split(':')[1]
                                f.source = ref_id
                                break
                    # replace SeqFeature with AutarkicSeqFeature in authoritative records list
                    strain_contig_records[isolate_id][record.id].features[i] = f
                    genes[isolate_id][ltag] = f
                    gene_names[isolate_id][ltag] = f.qualifiers['gene'][0] if 'gene' in f.qualifiers else ltag
                    unique_gene_names.add(gene_names[isolate_id][ltag])
                    node_data.append(prepare_node(isolate_id, ltag, f))


    iso = list(genes.keys())
    # Create log directories
    for sample_name in iso:
        os.makedirs(os.path.join(logs_dir, sample_name), exist_ok=True)

    for subdir in [logs_dir, reports_dir]:
        os.makedirs(subdir, exist_ok=True)

    # Initialize graphs used for collection correction data
    renames  = nx.Graph()
    renames.add_nodes_from(node_data)
    # https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.copy.html
    # TODO: We need the node attributes (gene names) to update in the new graph,
    # so this type of copy might not be appropriate. Double check.
    addition_references = renames.copy()
    additions_by_sample = defaultdict(set)

    # erase mcps file if it exists.
    # we're going to loop-write in append mode and we don't want to keep accumulating output from previous runs.
    if os.path.isfile(mcps_fn):
        open(mcps_fn, 'w').close()

    iso_pairs = list(itertools.combinations(range(len(iso)), r=2))

    with Pool(args.nproc) as p:
        mcp_results = p.map(mincanpairs, iso_pairs)

    # Create symlinks to existing log files
    for sample_name in iso:
        for other_sample in iso:
            if other_sample <= sample_name:  # Skip self and already handled pairs
                continue
            target_file = os.path.join(logs_dir, sample_name, f"{other_sample}.log")
            symlink_file = os.path.join(logs_dir, other_sample, f"{sample_name}.log")
            if os.path.exists(target_file) and not os.path.exists(symlink_file):
                rel_path = os.path.relpath(target_file, os.path.dirname(symlink_file))
                os.symlink(rel_path, symlink_file)

    for mcps, equivs, adds in mcp_results:
        with open(mcps_fn, 'a') as mcps_fh:
            for gene_pair in mcps:
                iso1, iso2 = mcps[gene_pair].keys()
                mcp1 = ' '.join(f"{sign(genes[iso1][g].location.strand)}{gene_names[iso1][g]}" for g in mcps[gene_pair][iso1]['segment'])
                mcp2 = ' '.join(f"{sign(genes[iso2][g].location.strand)}{gene_names[iso2][g]}" for g in mcps[gene_pair][iso2]['segment'])
                print('\t'.join([
                    gene_pair[0],
                    gene_pair[1],
                    iso1,
                    mcp1,
                    iso2,
                    mcp2,
                ]), file=mcps_fh)

                for entry in adds[gene_pair]:
                    # Make a new frozendict without the `ref` key, which will prevent
                    # duplicate additions from collapsing since the `ref` will never
                    # be the same
                    standard_entry = entry.delete('ref')
                    # Add (nonredundantly by coordinates) to our set of candidate additions for this strain
                    additions_by_sample[entry['sample']].add(standard_entry)
                    # Keep track of the strain reference for the addition in case its own name changes.
                    addition_references.add_edge(standard_entry, entry['ref'])

                renames.add_edges_from(equivs[gene_pair])

    # ensure transitivity of feature additions names.
    for sample in additions_by_sample:
        for addition in additions_by_sample[sample]:
            for node_id_pair in itertools.combinations(addition_references.neighbors(addition), 2):
                renames.add_edge(*node_id_pair)

    # Apply same criteria to resolving candidate renames as used for MCL postprocessing
    # TODO: we're working with the default orf prefix only currently.
    new_name_counter = designator.find_next_increment(unique_gene_names)
    name_changes = resolve_clusters(
        renames,
        new_name_counter,
        logfile=renames_fn,
    )

    additions_fh = open(additions_fn, 'w')
    for sample in additions_by_sample:
        sample_candidate_additions = postprocess_additions(
            additions_by_sample[sample],
            addition_references,
            strain_contig_records,
        )
        for contig in strain_contig_records[sample]:
            # log candidate additions pre-merge
            for feature in sample_candidate_additions[contig]:
                additions_fh.write("\t".join([
                    contig,
                    str(feature.location.start),
                    str(feature.location.end),
                    extractor.get_gene(feature),
                    dummy_bed_score,
                    sign(feature.location.strand),
                ]) + '\n')

            (_, _, _, _, _, _, G_overlap) = compare(
                sample_candidate_additions[contig],
                strain_contig_records[sample][contig].features,
                eliminate_colocated=False,
            )
            strain_contig_records[sample][contig].features = merge(G_overlap)

        if input_type != "bed":
            designator.assign_locus_tags(
                {
                    contig: strain_contig_records[sample][contig].features
                    for contig in strain_contig_records[sample]
                },
                prefix=sample,
            )
            out_gbk = os.path.join(args.outdir, f"{sample}.gbk")
            SeqIO.write(
                strain_contig_records[sample].values(),
                out_gbk,
                "genbank",
            )
            converter.convert_gbk_to_gff(out_gbk)
    additions_fh.close()

    if input_type == "bed":
        with open(os.path.join(args.outdir, "blocks_coords.bed"), 'w') as out_bed:
            for sample in strain_contig_records:
                for contig in strain_contig_records[sample]:
                    for feature in strain_contig_records[sample][contig].features:
                        out_bed.write("\t".join([
                            contig,
                            str(feature.location.start),
                            str(feature.location.end),
                            extractor.get_gene(feature),
                            dummy_bed_score,
                            sign(feature.location.strand),
                        ]) + '\n')


def overlap(loc1, loc2):
    loc1_start, loc1_end = loc1
    loc2_start, loc2_end = loc2
    overlap = (min(loc1_end, loc2_end) - max(loc1_start, loc2_start))
    if overlap > 0:
        return True
    else:
        return False

blast_fields = [
    'qseqid',
    'sseqid',
    'length',
    'pident',
    'qcovhsp',
    'evalue',
    'bitscore',
    'qstart',
    'qend',
    'sstart',
    'send',
    'qlen',
    'slen',
    'qseq',
    'sseq',
]
blast_outfmt   = '6 ' + ' '.join(blast_fields)

def summarize(blast_output):
    results = []
    blast_output = blast_output.strip()
    if not blast_output:
        return results
    for line in blast_output.strip().split('\n'):
        data = line.strip().split('\t')
        line_dict = dict(zip(blast_fields, data))
        for key in line_dict:
            if key in [
                    'qframe',
                    'sframe',
                    'qlen',
                    'slen',
            ]:
                line_dict[key] = int(line_dict[key])
            elif key in [
                    'pident',
            ]:
                line_dict[key] = float(line_dict[key])
        line_dict['qcov'] = (len(line_dict['qseq'].replace("-","")) / line_dict['qlen']) * 100
        line_dict['scov'] = (len(line_dict['sseq'].replace("-","")) / line_dict['slen']) * 100

        results.append(line_dict)
    return results

def tblastn_seg(qry_seg_genes, sub_loc, qry_ind, sub_ind):
    cmd = [
        'tblastn',
        '-db_gencode', str(cnf.genetic_code),
        '-outfmt',
        blast_outfmt,
        '-query', qry_seg_genes,
        '-subject', seq_files[iso[sub_ind]],
    ]
    ps = subprocess.run(cmd, check=True, capture_output=True, text=True)
    results = summarize(ps.stdout)
    return [
        r for r in results
        # see whether the coordinates of the best hit overlaps with subject segment
        # (we don't just blast against the segment for the case where any of the segment's genes overlap the flanking genes and thus go beyond the segment coordinates)
        if (
                overlap(sub_loc, (int(r['sstart']), int(r['send'])))
                and (
                    r['pident'] >= cnf.blast.min_identity
                    and (
                        r['scov'] >= cnf.blast.min_coverage
                        or r['qcov'] >= cnf.blast.min_coverage
                    )
                )
        )
    ]

def write_log_stanza(f, strain_pair, gene_pair, equivalences, seg1_additions, seg2_additions):
    """
    Log incremental corrections (from a single pairwise comparison) to a given log file.
    :param f: file handle of log file
    :param strain_pair: 2-tuple of strain IDs that were compared.
    :param gene_pair: 2-tuple of gene names of the minimum candidate pair
    :param equivalences: list of 2-tuples of locus tags that were marked as orthologs
    :param seg1_additions:
       list of dictionaries (generated by compare_segments) describing features to be added to strain 1 (from strain_pair).
    :param seg2_additions:
       as above, but for strain 2.
    """
    iso1, iso2 = strain_pair
    other_strain = {iso1: iso2, iso2:iso1}
    print(f"{gene_pair[0]}\t{gene_pair[1]}", file=f)
    for addn in seg1_additions + seg2_additions:
        gene_added = gene_names[other_strain[addn['sample']]][addn['ref']]
        print(f"//\t..\t..\tADD\t{addn['sample']}\t{gene_added}\t{addn['start']}\t{addn['end']}\t{addn['strand']}", file=f)
    for ltag1, ltag2 in equivalences:
        gene1 = gene_names[iso1][ltag1]
        gene2 = gene_names[iso2][ltag2]
        print(f"//\t..\t..\tRENAME\t{ltag1}|{gene1}\t{ltag2}|{gene2}", file=f)

def compare_segments(seg0_genes, seg1_genes, seg0_genes_file, seg1_genes_file, seg0_coords, seg1_coords, sample0_ind, sample1_ind):
    seg0_genes_unmatched = seg0_genes.copy()
    seg1_genes_unmatched = seg1_genes.copy()

    if not seg0_genes or not seg1_genes:
        bbh_results = {}
        seg0_qry_stats = []
        seg1_qry_stats = []
    else:
        bbh_results, seg0_qry_stats, seq1_qry_stats = BLAST.bidirectional_best_hit(
            seg0_genes_file,
            seg1_genes_file,
            identify=lambda seqid: seqid.split('|')[0],
            strict=False,
        )

    equivalences = []
    seg0_matches = {}
    seg1_matches = {}

    for ltag0 in bbh_results:
        if not bbh_results[ltag0]:
            continue

        gene0 = gene_names[iso[sample0_ind]][ltag0]
        ltag1 = bbh_results[ltag0]
        gene1 = gene_names[iso[sample1_ind]][ltag1]

        seg0_matches[ltag0] = True
        seg0_genes_unmatched.remove(ltag0)
        seg1_matches[ltag1] = True
        if ltag1 in seg1_genes_unmatched:
            seg1_genes_unmatched.remove(ltag1)

        # nothing to see here
        if gene0 == gene1:
            continue

        # can make this a 3-tuple with the third element being edge data
        # where we can keep blast results if we wanted to.
        equivalences.append((ltag0, ltag1))

    seg1_additions = []
    if seg0_genes_unmatched:
        seg0_genes_vs_intergene1 = tblastn_seg(seg0_genes_file, seg1_coords, sample0_ind, sample1_ind)
        for resultset in seg0_genes_vs_intergene1:
            ltag0, gene0, strand_of_0 = resultset['qseqid'].split('|')
            strand_of_0 = int(strand_of_0)
            if ltag0 in seg0_matches:
                continue
            seg0_matches[ltag0] = True
            start, end = int(resultset['sstart']), int(resultset['send'])
            # Check the strand of the hit and account for the stop codon position
            # since the protein sequence query doesn't represent it.
            # TODO: We should actually make sure this is a stop codon and not a run-on frame.
            #       will need to use stopseeker() for that.
            if start < end:
                sign_flip = 1
                end += 3
            else:
                sign_flip = -1
                start -= 3
            seg1_additions.append(frozendict({
                'sample':iso[sample1_ind],
                'contig_id':resultset['sseqid'],
                'start':min(start,end) - 1,
                'end':max(start,end),
                'strand':sign_flip,
                'ref':ltag0,
            }))
            seg0_genes_unmatched.remove(ltag0)

    seg0_additions = []
    if seg1_genes_unmatched:
        seg1_genes_vs_intergene0 = tblastn_seg(seg1_genes_file, seg0_coords, sample1_ind, sample0_ind)
        for resultset in seg1_genes_vs_intergene0:
            ltag1, gene1, strand_of_1 = resultset['qseqid'].split('|')
            strand_of_1 = int(strand_of_1)
            if ltag1 in seg1_matches:
                continue
            seg1_matches[ltag1] = True
            start, end = int(resultset['sstart']), int(resultset['send'])
            # Check the strand of the hit and account for the stop codon position
            # since the protein sequence query doesn't represent it.
            # TODO: We should actually make sure this is a stop codon and not a run-on frame.
            #       will need to use stopseeker() for that.
            if start < end:
                sign_flip = 1
                end += 3
            else:
                sign_flip = -1
                start -= 3
            seg0_additions.append(frozendict({
                'sample':iso[sample0_ind],
                'contig_id':resultset['sseqid'],
                'start':min(start,end) - 1,
                'end':max(start,end),
                'strand':sign_flip,
                'ref':ltag1,
            }))
            seg1_genes_unmatched.remove(ltag1)

    return equivalences, seg0_additions, seg1_additions, seg0_genes_unmatched, seg1_genes_unmatched

def get_segment(sample, pair):
    '''
    Get all genes occurring between the pair defining the interval.
    Note - this relies on all pair members being unique in their respective
           genome, but we have ensured that while creating the full list of pairs
    '''
    ltag_list = list(genes[sample].keys())
    coordA = ltag_list.index(pair[0])
    coordB = ltag_list.index(pair[1])
    return ltag_list[coordA:coordB+1]

def postprocess_additions(strain_additions, addition_refs, strain_contig_records):
    """
    Apply reference-based coordinate correction and resolve conflicting entries among the candidate additions for a strain.

    :param strain_additions: set of frozendicts representing coordinates to be added
    :param addition_refs: networkx Graph with frozendict nodes linked to strain reference locus tags with attached SeqFeature attributes (from prepare_node())
    :return: list of SeqFeature objects representing the postprocessed candidate additions that will need to be merged.
    """

    candidate_strain_additions = defaultdict(list)

    # retrieve the final name for the annotated sources
    for i, candidate in enumerate(strain_additions):
        candidate_id = f"L2_{i:05d}"
        candidate_feature = AutarkicSeqFeature(
            type='CDS',
            location=SimpleLocation(
                candidate['start'],
                candidate['end'],
                candidate['strand'],
                ref=candidate['contig_id'],
            ),
            references={
                candidate['contig_id']: strain_contig_records[candidate['sample']][candidate['contig_id']].seq,
            },
            qualifiers={
                'locus_tag':[candidate_id],
            }
        )
        candidate_feature.qualifiers['translation'] = [
            str(translate(candidate_feature.extract(), table=cnf.genetic_code))
        ]

        possible_names = defaultdict(set)
        for source_id in addition_refs.neighbors(candidate):
            source = addition_refs.nodes[source_id]
            possible_names[extractor.get_gene(source['annotation'], tryhard=False)].add(source_id)

        candidate_feature_personas = {}
        for name in possible_names:
            candidate_feature_personas[name] = deepcopy(candidate_feature)
            candidate_feature_personas[name].qualifiers['gene'] = [name]
            ref_instance = addition_refs.nodes[
                next(iter(possible_names[name]))
            ]['annotation']
            if ref_instance.source:
                ref_feature_origin = ref_instance.source
            else:
                ref_feature_origin = ref_instance.location.parts[0].ref
            if not designator.is_reference(name):
                annomerge.ref_annotation[
                    designator.key_ref_gene(ref_feature_origin, name)
                ] = ref_instance
            ref_feature = ref_annotation[designator.key_ref_gene(ref_feature_origin, name)]
            candidate_feature_personas[name].source = ref_feature_origin

            # coordinate correction
            pseudoscan(
                candidate_feature_personas[name],
                ref_feature,
                attempt_rescue=True,
            )

        remaining_possible_names = list(possible_names.keys())
        while len(remaining_possible_names) > 1:
            # thunderdome tournament. Only one can prevail.
            name_contender1, name_contender2 = remaining_possible_names[0:2]
            include1, include2, evid, remark = thunderdome(
                candidate_feature_personas[name_contender1],
                candidate_feature_personas[name_contender2],
            )
            if not include1:
                remaining_possible_names.remove(name_contender1)
            else:
                remaining_possible_names.remove(name_contender2)
        else:
            candidate_strain_additions[candidate['contig_id']].append(
                candidate_feature_personas[remaining_possible_names[0]]
            )

    for contig in candidate_strain_additions:
        (_, _, _, G_overlaps) = cross_examine(sort_features(candidate_strain_additions[contig]))
        candidate_strain_additions[contig] = merge(G_overlaps)

    return candidate_strain_additions

def mincanpairs(sample_ind_pair):
    sample0_ind, sample1_ind = sample_ind_pair
    sample0, sample1 = iso[sample0_ind], iso[sample1_ind]

    # Determine log file path (always use lexicographically first sample as directory)
    if sample0 < sample1:
        log_file = os.path.join(logs_dir, sample0, f"{sample1}.log")
    else:
        log_file = os.path.join(logs_dir, sample1, f"{sample0}.log")
    log_fh = open(log_file, 'w')


    can_pairs = []
    i = 0

    iso0_ltags = list(genes[iso[sample0_ind]].keys())
    iso0_genes = list(gene_names[iso[sample0_ind]].values())
    iso1_ltags = list(genes[iso[sample1_ind]].keys())
    iso1_genes = list(gene_names[iso[sample1_ind]].values())

    # find all candidate pairs
    for ltag1a in iso0_ltags:
        gene1 = gene_names[iso[sample0_ind]][ltag1a]
        # stop when we get 1 gene from the end since there won't be a pair for it
        # (subtract 2 instead of 1 because of 0-indexing)
        if i == len(iso0_genes)-2:
            break
        # treat any gene with multiple copies as potentially different
        # because it will be complicated to track down which copy corresponds
        # to which between the two isolates.
        if iso0_genes.count(gene1)>1 or \
           gene1 not in iso1_genes or \
           iso1_genes.count(gene1)>1:
            i += 1
            continue
        ltag2a_pos = iso1_genes.index(gene1)
        ltag2a = iso1_ltags[ltag2a_pos]
        for ltag1b in iso0_ltags[i+2:]:
            gene2 = gene_names[iso[sample0_ind]][ltag1b]
            if iso0_genes.count(gene2)>1 or \
               gene2 not in iso1_genes or \
               iso1_genes.count(gene2)>1:
                continue
            ltag2b_pos = iso1_genes.index(gene2)
            ltag2b = iso1_ltags[ltag2b_pos]
            if ltag2a_pos > ltag2b_pos:
                (ltag2a, ltag2b) = (ltag2b, ltag2a)
            can_pairs.append([(gene1, gene2), (ltag1a,ltag1b), (ltag2a, ltag2b)])
            break
        i += 1

    # prepare final results by checking minimality
    mcps = {}
    equivs = {}
    additions = {}
    for gene_pair, iso0_pair, iso1_pair in can_pairs:
        segment1 = get_segment(iso[sample0_ind], iso0_pair)
        segment1_gene_names = [gene_names[iso[sample0_ind]][g] for g in segment1]
        segment2 = get_segment(iso[sample1_ind], iso1_pair)
        segment2_gene_names = [gene_names[iso[sample1_ind]][g] for g in segment2]
        # There's simply a deletion
        # TODO -- get rid of this -- could be unannotated in between!
        #if len(segment1)==2 or len(segment2)==2:
        #    continue
        if (
                (
                    segment1_gene_names[0] == segment2_gene_names[0]
                    and segment1_gene_names[-1] == segment2_gene_names[-1]
                    and segment1_gene_names[0:2] != segment2_gene_names[0:2]
                    and segment1_gene_names[-2:] != segment2_gene_names[-2:]
                ) or ( #inverted pair
                    segment1_gene_names[0] == segment2_gene_names[-1]
                    and segment1_gene_names[-1] == segment2_gene_names[0]
                    and segment1_gene_names[0:2] != segment2_gene_names[-2:][::-1]
                    and segment1_gene_names[-2:] != segment2_gene_names[0:2][::-1]
                )
        ):
            seg1_coords = (
                genes[iso[sample0_ind]][iso0_pair[0]].location.end,
                genes[iso[sample0_ind]][iso0_pair[1]].location.start + 1
            )
            seg2_coords = (
                genes[iso[sample1_ind]][iso1_pair[0]].location.end,
                genes[iso[sample1_ind]][iso1_pair[1]].location.start + 1
            )
            with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=cnf.tmpdir) as seg1_geneseqs, tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=cnf.tmpdir) as seg2_geneseqs:
                seg1_gene_records = [gene_seqs[iso[sample0_ind]][g] for g in segment1]
                SeqIO.write(seg1_gene_records, seg1_geneseqs, "fasta")
                seg2_gene_records = [gene_seqs[iso[sample1_ind]][g] for g in segment2]
                SeqIO.write(seg2_gene_records, seg2_geneseqs, "fasta")
                (
                    equivalences,
                    seg1_additions,
                    seg2_additions,
                    seg1_unmatched,
                    seg2_unmatched,
                ) = compare_segments(
                    segment1,
                    segment2,
                    seg1_geneseqs.name,
                    seg2_geneseqs.name,
                    seg1_coords,
                    seg2_coords,
                    sample0_ind,
                    sample1_ind,
                )
            mcps[gene_pair] = {
                iso[sample0_ind]: {
                    'segment': segment1
                },
                iso[sample1_ind]: {
                    'segment': segment2
                }
            }
            additions[gene_pair] = seg1_additions + seg2_additions
            equivs[gene_pair] = equivalences

            write_log_stanza(
                log_fh,
                strain_pair=(iso[sample0_ind], iso[sample1_ind]),
                gene_pair=gene_pair,
                equivalences=equivalences,
                seg1_additions=seg1_additions,
                seg2_additions=seg2_additions,
            )

    log_fh.close()

    return mcps, equivs, additions
