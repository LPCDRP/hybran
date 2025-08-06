from copy import deepcopy
import os

from .bio import (
    sort_features,
    SeqIO,
)
from .lumberjack import load_fusions
from . import (
    converter,
    fileManager,
)


def main(args):
    annotations = fileManager.file_list(
        [args.annotations_dir],
        file_type='genbank',
    )

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except:
            sys.exit("Could not create output directory " + args.output)

    for annotation in annotations:
        basefilename = os.path.basename(annotation)
        name = os.path.splitext(basefilename)[0]
        fusion_report = os.path.join(args.annotations_dir, name, 'fusion_report.tsv')
        records = []
        for record in SeqIO.parse(annotation, "genbank"):
            fusions = load_fusions(fusion_report)
            record.features = defuse(record.features, fusions)
            records.append(record)
        out_annotation = os.path.join(args.output, basefilename)
        SeqIO.write(records, out_annotation, "genbank")
        converter.convert_gbk_to_gff(out_annotation)

def defuse(feature_list, fusions_by_ltag):
    out_features = []

    for feature in feature_list:
        if feature.type == 'gene':
            continue
        if (
                'locus_tag' in feature.qualifiers
                and feature.qualifiers['locus_tag'][0] in fusions_by_ltag
        ):
            ltag = feature.qualifiers['locus_tag'][0]
            fusion_type = fusions_by_ltag[ltag]['type']
            components = fusions_by_ltag[ltag]['components']
            subltags = (fusion_type == 'partial' or (fusion_type == 'whole' and len(components) > 2))
            sub = 0
            for component in components:
                if not 'locus_tag' in component.qualifiers:
                    if subltags:
                        component.qualifiers['locus_tag'] = [ f"{ltag}{chr(ord('A')+sub)}" ]
                        sub += 1
                    else:
                        component.qualifiers['locus_tag'] = [ltag]
                out_features.append(component)
        else:
            out_features.append(feature)

    out_features = sort_features(out_features)

    return out_features
