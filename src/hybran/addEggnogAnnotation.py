
import os
from Bio import SeqIO
from . import converter
from . import config
from . import designator

def parse_eggnog(ref_tax_ids):
    hybran_tmp_dir = config.hybran_tmp_dir
    orthologs_annotation = {}
    annotation_dict = {}
    orthologs_dict = {}
    hmm_generics = []

    for record in SeqIO.parse(hybran_tmp_dir + '/eggnog_seqs.fasta', 'fasta'):
        hmm_generics.append(record.id)
    with open(hybran_tmp_dir + '/eggnog-mapper-annotations/generics_diamond.emapper'
                              '.annotations', 'r') as diamond_results_annotation:
        for line in diamond_results_annotation:
            if line.startswith('#'):
                continue
            line_elements = line.strip().split('\t')
            if line_elements[0] in hmm_generics:
                continue
            if len(line_elements) < 13:
                annotation = 'NA'
            else:
                annotation = line_elements[10]
            if not any([_ in line_elements[1] for _ in ref_tax_ids]):
                orthologs_annotation[line_elements[0]] = annotation
                continue
            dict_val = ('Eggnog:diamond-seed', line_elements[1], annotation)
            annotation_dict[line_elements[0]] = dict_val

    annotated_generics = list(annotation_dict.keys())
    with open(hybran_tmp_dir +
              '/eggnog-mapper-annotations/generics_diamond.emapper.orthologs', 'r') as diamond_results_orthologs:
        for line in diamond_results_orthologs:
            line_elements = line.strip().split('\t')
            if line_elements[0] in annotated_generics:
                continue
            if len(line_elements) < 2:
                orthologs_dict[line_elements[0]] = ''
            else:
                orthologs = line_elements[1].split(',')
                orthologs_dict[line_elements[0]] = orthologs
                corresponding_ref = ''
                for ort in orthologs:
                    if any([_ in ort for _ in ref_tax_ids]):
                        corresponding_ref = ort
                        break
                if len(corresponding_ref) > 0:
                    annotation_dict[line_elements[0]] = ('Eggnog:diamond-ortholog',
                                                         corresponding_ref,
                                                         orthologs_annotation[line_elements[0]])
    with open(hybran_tmp_dir + '/eggnog-mapper-annotations/generics_hmm.emapper.annotations',
              'r') as hmm_results_annotation:
        for line in hmm_results_annotation:
            if line.startswith('#'):
                continue
            line_elements = line.strip().split('\t')
            if len(line_elements) < 13:
                annotation = 'NA'
            else:
                annotation = line_elements[10]
            if not any([_ in line_elements[1] for _ in ref_tax_ids]):
                orthologs_annotation[line_elements[0]] = annotation
                continue
            dict_val = ('Eggnog:hmm-seed', line_elements[1], annotation)
            annotation_dict[line_elements[0]] = dict_val

    annotated_generics = list(annotation_dict.keys())
    with open(hybran_tmp_dir + '/eggnog-mapper-annotations/generics_hmm.emapper.orthologs', 'r') as hmm_results_orthologs:
        for line in hmm_results_orthologs:
            line_elements = line.strip().split('\t')
            if line_elements[0] in annotated_generics:
                continue
            if len(line_elements) < 2:
                orthologs_dict[line_elements[0]] = ''
            else:
                orthologs = line_elements[1].split(',')
                orthologs_dict[line_elements[0]] = orthologs
                corresponding_ref = ''
                for ort in orthologs:
                    if any([_ in ort for _ in ref_tax_ids]):
                        corresponding_ref = ort
                        break
                if len(corresponding_ref) > 0:
                    annotation_dict[line_elements[0]] = ('Eggnog:hmm-ortholog',
                                                         corresponding_ref,
                                                         orthologs_annotation[line_elements[0]])
    return annotation_dict


def update_gbks(ref_tax_ids, ref_gene_dict):
    """
    Update isolate annotation files with additional reference gene names
    mapped via eggnog

    :param ref_tax_ids: list of NCBI taxonomy IDs of the references
    :param ref_gene_dict: dict mapping reference locus tags to gene names
    """

    annotations_to_add = parse_eggnog(ref_tax_ids=ref_tax_ids)
    gbks = [gbk for gbk in os.listdir(os.getcwd()) if gbk.endswith('.gbk')]
    for gbk in gbks:
        gbk_fp = gbk
        output_fp = gbk
        output_recs = []
        reference_generics_dict = {}
        generic_genes_in_isolate = {}
        for record in SeqIO.parse(gbk_fp, 'genbank'):
            rec_features = record.features
            for feature in rec_features:
                if feature.type != 'CDS':
                    continue
                if not designator.is_unannotated(feature.qualifiers['gene'][0]):
                    continue
                if feature.qualifiers['gene'][0] not in annotations_to_add.keys():
                    continue
                annotation_info = annotations_to_add[feature.qualifiers['gene'][0]]
                ref_gene_id = annotation_info[1].split('.')[1]
                generic_genes_in_isolate[feature.qualifiers['gene'][0]] = ref_gene_id
                if ref_gene_id in reference_generics_dict.keys():
                    reference_generics_dict[ref_gene_id].append(feature.qualifiers['gene'][0])
                else:
                    reference_generics_dict[ref_gene_id] = [feature.qualifiers['gene'][0]]
                new_note = annotation_info[0] + ':' + annotation_info[1].split('.')[1]
                if 'note' not in feature.qualifiers.keys():
                    feature.qualifiers['note'] = [new_note]
                else:
                    feature.qualifiers['note'].append(new_note)
                # If it was transferred from a genome that was already annotated with EggNOG
                # Skip this CDS
                if [n for n in feature.qualifiers['note'] if n.startswith('Eggnog')]:
                    continue
                if annotation_info[2] != 'NA':
                    annotation_note = 'Eggnog:annotation:' + annotation_info[2]
                    feature.qualifiers['note'].append(annotation_note)
            for feature in rec_features:
                if feature.type != 'CDS':
                    continue
                if feature.qualifiers['gene'][0] not in generic_genes_in_isolate.keys():
                    continue
                else:
                    annotation_info = annotations_to_add[feature.qualifiers['gene'][0]]
                    eggnog_gene = annotation_info[1].split('.')[1]
                    eggnog_annotation = annotation_info[2]
                    corresponding_ref = generic_genes_in_isolate[feature.qualifiers['gene'][0]]
                    num_of_unique_generics = set(reference_generics_dict[corresponding_ref])
                    if len(num_of_unique_generics) > 1:
                        feature.qualifiers['note'].append(
                            'This reference gene is possibly split across multiple CDSs: Multiple generic ORFs '
                            'in this isolate annotated with same reference gene ID by EggNOG')
                    else:
                        try:
                            feature.qualifiers['gene'][0] = ref_gene_dict[eggnog_gene]
                        except KeyError:
                            feature.qualifiers['gene'][0] = eggnog_gene
            output_recs.append(record)
        SeqIO.write(list(output_recs), output_fp, 'genbank')
