#!/usr/bin/env python2.7

from Bio import SeqIO
from sets import Set
import pickle
import os


GROUPHOME = os.environ['GROUPHOME']


def main():
    mtb_pickle = pickle.load(open('isolate_gene_name.pickle', 'rb'))
    gbk_files = os.listdir(GROUPHOME + '/data/annotation/')
    isolates = [isolate_id.split('.')[0] for isolate_id in gbk_files if isolate_id.endswith('.gbk')]
    for isolate in mtb_pickle.keys():
        print(isolate)
        if isolate not in isolates:
            print('Isolate ' + isolate + ' absent in $GROUPHOME/data/annotation')
        else:
            genbank_file = GROUPHOME + '/data/annotation/' + isolate + '.gbk'
            isolate_records = list(SeqIO.parse(genbank_file, 'genbank'))
            update_mtb_dict = mtb_pickle[isolate]
            locus_to_update = update_mtb_dict.keys()
            modified_locus = []
            rec_num = 1
            for rec in isolate_records:
                if rec_num > 1:
                    break
                rec_num += 1
                for feature in rec.features:
                    if feature.type != 'CDS':
                        continue
                    elif feature.qualifiers['locus_tag'][0] in locus_to_update and feature.qualifiers['locus_tag'][0] not in modified_locus:
                        #print(feature.qualifiers['locus_tag'][0])
                        #print(update_mtb_dict[feature.qualifiers['locus_tag'][0]])
                        #print(feature.qualifiers['gene'])
                        if feature.qualifiers['gene'][0].startswith('L'):
                            feature.qualifiers['gene'][0] = update_mtb_dict[feature.qualifiers['locus_tag'][0]]
                        elif feature.qualifiers['gene'][0] == update_mtb_dict[feature.qualifiers['locus_tag'][0]]:
                            modified_locus.append(feature.qualifiers['locus_tag'][0])
                            continue
                        else:
                            print('Discordant assignment of gene name')
                            print('Original gene name: ' + feature.qualifiers['gene'][0])
                            print('New gene name: ' + update_mtb_dict[feature.qualifiers['locus_tag'][0]])
                            if 'gene_synonym' in feature.qualifiers.keys():
                                feature.qualifiers['gene_synonym'].append(update_mtb_dict[feature.qualifiers['locus_tag'][0]])
                            else:
                                feature.qualifiers['gene_synonym'] = [update_mtb_dict[feature.qualifiers['locus_tag'][0]]]
                        modified_locus.append(feature.qualifiers['locus_tag'][0])
                    elif feature.qualifiers['locus_tag'][0] in locus_to_update and feature.qualifiers['locus_tag'][0] in modified_locus:
                        print('ERROR: Updated locus tag previously')
                    else:
                        continue
                if len(Set(locus_to_update).intersection(Set(modified_locus))) < len(locus_to_update):
                    print('The following locus_tags are missing in the genbank file')
                    print(Set(locus_to_update).difference(Set(modified_locus)))
            SeqIO.write(isolate_records, genbank_file, 'genbank')


if __name__ == '__main__':
    main()


