
import os
from Bio import SeqIO
from . import converter
from . import config
from . import designator

def parse_eggnog(ref_tax_id):
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
                annotation = line_elements[-1]
            if ref_tax_id not in line_elements[1]:
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
                corresponding_rv = ''
                for ort in orthologs:
                    if ref_tax_id in ort:
                        corresponding_rv = ort
                        break
                if len(corresponding_rv) > 0:
                    annotation_dict[line_elements[0]] = ('Eggnog:diamond-ortholog',
                                                         corresponding_rv,
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
                annotation = line_elements[-1]
            if ref_tax_id not in line_elements[1]:
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
                corresponding_rv = ''
                for ort in orthologs:
                    if ref_tax_id in ort:
                        corresponding_rv = ort
                        break
                if len(corresponding_rv) > 0:
                    annotation_dict[line_elements[0]] = ('Eggnog:hmm-ortholog',
                                                         corresponding_rv,
                                                         orthologs_annotation[line_elements[0]])
    return annotation_dict


def update_gbks(script_dir, ref_tax_id):
    """
    Update isolate annotation files with additional reference gene names
    mapped via eggnog

    :param ref_tax_id: NCBI taxonomy ID of the reference
    """

    def h37rv_gene_names():
        pe_ppe_rv = []
        pe_ppe_gene = []
        genes_dict = {}
        rv_genes = os.sep.join([script_dir,'resources','tuberculist_genes.tsv'])
        with open(rv_genes, 'r') as gene_file:
            for line in gene_file:
                column = line.rstrip('\n').split('\t')
                if 'locus_tag' in column[0]:
                    header = column
                    continue
                genes_dict[column[0]] = column[1]
                if 'PE/PPE' in column[header.index('functional_category')]:
                    pe_ppe_rv.append(column[0])
                    pe_ppe_gene.append(column[1])
        return pe_ppe_rv, pe_ppe_gene, genes_dict
    annotations_to_add = parse_eggnog(ref_tax_id=ref_tax_id)
    pe_ppe_rv_list, pe_ppe_gene_list, h37rv_genes = h37rv_gene_names()
    gbks = [gbk for gbk in os.listdir(os.getcwd()) if gbk.endswith('.gbk')]
    for gbk in gbks:
        gbk_fp = gbk
        output_fp = gbk
        counter = 1
        output_recs = []
        reference_generics_dict = {}
        generic_genes_in_isolate = {}
        for record in SeqIO.parse(gbk_fp, 'genbank'):
            if counter > 1:
                output_recs.append(record)
                counter += 1
                continue
            rec_features = record.features
            for feature in rec_features:
                if feature.type != 'CDS':
                    continue
                if not designator.is_unannotated(feature.qualifiers['gene'][0]):
                    continue
                if feature.qualifiers['gene'][0] not in annotations_to_add.keys():
                    continue
                annotation_info = annotations_to_add[feature.qualifiers['gene'][0]]
                rv_gene_id = annotation_info[1].split('.')[1]
                generic_genes_in_isolate[feature.qualifiers['gene'][0]] = rv_gene_id
                if rv_gene_id in reference_generics_dict.keys():
                    reference_generics_dict[rv_gene_id].append(feature.qualifiers['gene'][0])
                else:
                    reference_generics_dict[rv_gene_id] = [feature.qualifiers['gene'][0]]
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
                    corresponding_rv = generic_genes_in_isolate[feature.qualifiers['gene'][0]]
                    num_of_unique_generics = set(reference_generics_dict[corresponding_rv])
                    if len(num_of_unique_generics) > 1:
                        feature.qualifiers['note'].append(
                            'This reference gene is possibly split across multiple CDSs: Multiple generic ORFs '
                            'in this isolate annotated with same reference gene ID by EggNOG')
                    elif 'PE' not in eggnog_annotation and \
                            'PPE' not in eggnog_annotation and \
                            eggnog_gene not in pe_ppe_rv_list:
                        try:
                            feature.qualifiers['gene'][0] = h37rv_genes[eggnog_gene]
                        except KeyError:
                            feature.qualifiers['gene'][0] = eggnog_gene
            output_recs.append(record)
            counter += 1
        SeqIO.write(list(output_recs), output_fp, 'genbank')
