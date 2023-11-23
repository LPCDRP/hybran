import os

from collections import defaultdict

from . import fileManager
from . import annomerge
from . import designator
from .bio import SeqIO

def main(args):

    gbk_file1 = args.annotations[0]
    gbk_file2 = args.annotations[1]
    outfile = args.output

    results, pseudo_results = compare(gbk_file1, gbk_file2)

    nl = "\n"
    statement = (
        f"Overlaps: {len(results['overlap'])} | Exact: {len(results['exact'])} |" +
        f"X_unique: {len(results['x_unique'])} | Y_unique: {len(results['y_unique'])}{nl}" +
        f"Pseudo Overlaps: {len(pseudo_results['overlap'])} | Pseudo Exact: {len(pseudo_results['exact'])} |" +
        f"Pseudo_X_unique: {len(pseudo_results['x_unique'])} | Pseudo_Y_unique: {len(pseudo_results['y_unique'])}"
    )
    print(statement)
    breakpoint()

def furthest_location(loc1, loc2):
    """
    Returns the SeqFeature.Location that is farthest downstream or returns None if they are exact matches.
    :param loc1: SeqFeature location.
    :param loc2: SeqFeature location.
    """
    same_start = False
    same_stop = False
    if loc1.start == loc2.start:
        same_start = True
    if loc1.end == loc2.end:
        same_stop = True

    if same_start and same_stop:
        return None
    elif not same_start and same_stop:
        if loc1.start > loc2.start:
            return loc1
        return loc2
    else:
        if loc1.end > loc2.end:
            return loc1
        return loc2

def produce_record(gbk):
    """
    Populate a defaultdictionary with all of the features from an annotation file

    :param gbk: String of a path to a .gbk annotation file.
    :return anno_dict: Defaultdictionary with all CDS features from the input .gbk file.
    :return anno_list: List of gene names from the input.gbk file ordered by position.
    """
    anno_dict = defaultdict(list)
    anno_list = []
    anno_record = SeqIO.read(gbk, "genbank")
    for f in anno_record.features:
        if f.type == 'CDS' or f.type == 'pseudo':
            gene_name = (f.qualifiers.get('locus_tag')[0] if f.qualifiers.get('gene') == None else f.qualifiers.get('gene')[0])
            if gene_name == 'mamB' and "RATT" in f.qualifiers['inference'][0]:
                continue
            anno_dict[gene_name].append(f)
            anno_list.append(gene_name)

    return anno_dict, anno_list
    

def compare(gbk1, gbk2):
    """
    Compare between two .gbk annotation files and produce a summary report
    
    :param gbk1: String of a path to a .gbk annotation file
    :param gbk2: String of a path to a .gbk annotation file
    :return results: Nested dictionary of the comparison results
    :return pseudo_results: Nested dictionary of the comparison results of only pseudo genes
    """

    results = defaultdict(list)
    pseudo_results = defaultdict(list)
    
    x_anno_dict, x_anno_list = produce_record(gbk1)
    y_anno_dict, y_anno_list = produce_record(gbk2)

    x_dupes = []
    y_dupes = []

    x = 0
    y = 0

    for i in range(x, len(x_anno_list)):
        x_gene = x_anno_list[i]
        x_dupes.append(x_gene)
        #Keeps track of which duplicate gene is being referenced.
        x_index = x_dupes.count(x_gene) - 1
        x_feature = x_anno_dict[x_gene][x_index]
        
        for j in range(y, len(y_anno_list)):
            y_gene = y_anno_list[j]
            y_dupes.append(y_gene)
            y_index = y_dupes.count(y_gene) - 1
            y_feature = y_anno_dict[y_gene][y_index]

            farthest_feature = furthest_location(x_feature.location, y_feature.location)

            #x_feature is ahead or both are the in the exact same location
            if farthest_feature == x_feature.location or farthest_feature == None:

                #exact same
                if farthest_feature == None:
                    results['exact'].append([x_feature, y_feature])
                    if designator.is_pseudo(x_feature.qualifiers) or designator.is_pseudo(y_feature.qualifiers):
                        pseudo_results['exact'].append([x_feature, y_feature])
                    
                #not exact same but overlapping
                elif annomerge.overlap_inframe(x_feature.location, y_feature.location):
                     results['overlap'].append([x_feature, y_feature])
                     if designator.is_pseudo(x_feature.qualifiers) or designator.is_pseudo(y_feature.qualifiers):
                        pseudo_results['overlap'].append([x_feature, y_feature])
                     
                #not exact same, no overlap, and the x_feature is ahead
                else:
                    results['x_unique'].append(x_feature)
                    if designator.is_pseudo(x_feature.qualifiers):
                        pseudo_results['x_unique'].append(x_feature)
                        
                if y+1 == len(y_anno_list):
                    y = len(y_anno_list)-1
                else:
                    y = y + 1

                    #dont need to pop the feature from x_dupes because "continue" will only loop through the 'y' section
                    continue
                        
            #Y_feature is ahead
            elif farthest_feature == y_feature.location:

                #need to pop feature from y_dupes because "break" will start loop back in the 'x' section and we don't want
                #to append the y_feature twice.
                y_dupes.pop()

                #not exact same but overlapping
                if annomerge.overlap_inframe(x_feature.location, y_feature.location):
                    results['overlap'].append([x_feature, y_feature])
                    if designator.is_pseudo(x_feature.qualifiers) or designator.is_pseudo(y_feature.qualifiers):
                        pseudo_results['overlap'].append([x_feature, y_feature])
                    
                #not exact same, no overlap, y_feature is ahead
                else:
                    results['y_unique'].append([x_feature, y_feature])
                    if designator.is_pseudo(y_feature.qualifiers):
                        pseudo_results['y_unique'].append(x_feature)
                    
                if x + 1 == len(x_anno_list):
                    x = len(x_anno_list)-1
                else:
                    x = x + 1
                break
                
            else:
                print("something is wrong, shouldn't make it here.")

    return results, pseudo_results
