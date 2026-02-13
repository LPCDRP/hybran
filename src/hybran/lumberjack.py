import Bio
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord

from .bio import AutarkicSeqFeature, SeqIO
from . import designator
from .extractor import get_ltag, get_gene
from .util import nacast


def log_cluster_resolution(res_data, logfile):
    header = [
        'cluster_id',
        'cluster_type',
        'cluster_size',
        'locus_tag',
        'original_gene_name',
        'final_gene_name',
    ]
    print('\t'.join(header), file=logfile)
    none_last = lambda x: x if x is not None else 'zzz'
    for resolution in sorted(res_data, key=lambda x: (x[0], none_last(x[1]), nacast(x[5]))):
        print(
            '\t'.join([nacast(field) for field in resolution]),
            file=logfile,
        )

def log_feature_fate(feature, logfile, superior=None, evid=None, remark=None):
    """
    General-purpose logging function to print out a gene's information and a comment
    :param feature: A SeqFeature object
    :param logfile: An open filehandle
    :param superior: SeqFeature representing the superior annotation displacing feature.
    :param evid: str evidence code
    :param remark: str explanation
    """
    if superior is not None:
        superior_locus_tag = get_ltag(superior)
        superior_gene_name = get_gene(superior)
    else:
        superior_locus_tag = superior_gene_name = None

    locus_tag = get_ltag(feature)
    gene_name = get_gene(feature)

    row = [
        locus_tag,
        gene_name,
        superior_locus_tag,
        superior_gene_name,
        evid,
        remark,
    ]
    print(
        '\t'.join([nacast(field) for field in row]),
        file=logfile
    )

def log_feature_fates(info_list, logfile):
    """
    Calls log_feature_fate on each element in info_list and creates headers
    :param info_list:
      list of dictionaries with required key/values
        - 'feature': SeqFeature that is rejected
        - 'evid': str evidence code
      and optional key/values
        - 'remark': str explanation
        - 'superior': SeqFeature of prevailing annotation
    :param logfile: An open filehandle
    """
    header = [
        'locus_tag',
        'gene_name',
        'rival_locus_tag',
        'rival_gene_name',
        'evidence_codes',
        'remark',
    ]
    print('\t'.join(header), file=logfile)
    if info_list:
        info_list.sort(key=lambda _:_['feature'].qualifiers['locus_tag'][0])

    [log_feature_fate(**_, logfile=logfile) for _ in info_list]

def log_coord_correction(feature, logfile):
    """
    This function is used to log gene information when coord_check() fixes a start/stop postition for individual features.
    It is called by log_coord_corrections().
    :param feature: AutarkicSeqFeature object
    :param logfile: An open filehandle
    """
    locus_tag = feature.qualifiers['locus_tag'][0]
    gene_name = feature.qualifiers['gene'][0]
    strand = str(feature.location.strand)
    og_start = (int(feature.og.location.start) + 1)
    og_end = (int(feature.og.location.end))
    new_start = (int(feature.corr.location.start) + 1)
    new_end = (int(feature.corr.location.end))
    start_fixed = str(og_start != new_start).lower()
    stop_fixed = str(og_end != new_end).lower()
    gene_length_ratio = f"{(og_end - og_start)/(new_end - new_start):.3f}"

    if feature.location.strand == -1:
        start_fixed, stop_fixed = stop_fixed, start_fixed

    line = [
        locus_tag,
        gene_name,
        strand,
        og_start,
        og_end,
        new_start,
        new_end,
        start_fixed,
        stop_fixed,
        gene_length_ratio,
        feature.corr_accepted,
    ]
    print('\t'.join(str(v) for v in line), file=logfile)

def log_coord_corrections(features_by_contig_dict, logfile):
    """
    Log status of all correctable features.
    :param features_by_contig_dict: dict of lists of AutarkicSeqFeatures where key is contig name
    :param logfile: an open filehandle
    """
    header = [
        'locus_tag',
        'gene_name',
        'strand',
        'og_start',
        'og_end',
        'new_start',
        'new_end',
        'fixed_start_codon',
        'fixed_stop_codon',
        'gene_length_ratio',
        'accepted',
    ]
    print('\t'.join(header), file=logfile)

    for contig in features_by_contig_dict:
        for feature in features_by_contig_dict[contig]:
            if feature.corr_possible:
                log_coord_correction(feature, logfile)

def log_pseudos(features_by_contig_dict, logfile):
    """
    Function to log the pseudo status of all anomalous features determined by pseudoscan()
    :param features_by_contig_dict: dict of lists of AutarkicSeqFeatures where key is contig name
    :param logfile: an open filehandle
    """
    #Note codes are categorized by a '0' or '1' and correspond to 'False' and 'True' respectively
    #D3 = Divisible by three [0/1]
    #VS = Valid start [0/1]
    #VE = Valid end [0/1]
    #RCS = Reference corresponding start [0/1]
    #RCE = Reference corresponding end [0/1]
    #BOK = Blast OK [0/1]

    header = [
        'locus_tag',
        'gene_name',
        'pseudo',
        'evidence_codes',
        'summmary',
        'div_by_3',
        'valid_start_codon',
        'valid_stop_codon',
        'ref_corr_start',
        'ref_corr_end',
        'blast_ok',
    ]
    print('\t'.join(header), file=logfile)
    for contig in features_by_contig_dict:
        for feature in features_by_contig_dict[contig]:
            if feature.ps_evid:
                summary = ';'.join([
                    f"D3{nacast(feature.d3)}",
                    f"VS{nacast(feature.vs)}",
                    f"VE{nacast(feature.ve)}",
                    f"RCS{nacast(feature.rcs)}",
                    f"RCE{nacast(feature.rce)}",
                    f"BOK{nacast(feature.bok)}",
                ])
                line = [
                    feature.qualifiers['locus_tag'][0],
                    feature.qualifiers['gene'][0],
                    nacast(designator.is_pseudo(feature.qualifiers)),
                    ';'.join(feature.ps_evid),
                    summary,
                    nacast(feature.d3),
                    nacast(feature.vs),
                    nacast(feature.ve),
                    nacast(feature.rcs),
                    nacast(feature.rce),
                    nacast(feature.bok),
                ]
                print('\t'.join(str(v) for v in line), file=logfile)

def log_fusions(features_by_contig_dict, logfile):
    header = [
        'fusion_type',
        'locus_tag',
        'gene_name',
        'strand',
        'start',
        'end',
    ]
    print('\t'.join(header), file=logfile)
    for contig in features_by_contig_dict:
        for feature in features_by_contig_dict[contig]:
            if feature.fusion_type:
                gene_info = []
                for f in [feature] + feature.fusion_components:
                    gene_info.append([
                        '..',
                        f.qualifiers['locus_tag'][0] if 'locus_tag' in f.qualifiers else '.',
                        get_gene(f),
                        str(f.location.strand),
                        str(f.location.start + 1),
                        str(f.location.end),
                    ])
                gene_info[0][0] = feature.fusion_type
                for line in gene_info:
                    print('\t'.join(line), file=logfile)


def load_fusions(fusion_report_file):
    fusions_by_ltag = {}
    with open(fusion_report_file, 'r') as fusion_report:
        next(fusion_report)
        for line in fusion_report:
            line = line.strip()
            if not line.startswith('..'):
                [
                    fusion_type,
                    ltag,
                    gene,
                    strand,
                    start,
                    end,
                ] = line.split('\t')
                fusions_by_ltag[ltag] = {
                    'type': fusion_type,
                    'gene': gene,
                    'strand': strand,
                    'start': start,
                    'end': end,
                    'components': []
                }
            else:
                [
                    _,
                    comp_ltag,
                    comp_gene,
                    comp_strand,
                    comp_start,
                    comp_end,
                ] = line.split('\t')
                comp = Bio.SeqFeature.SeqFeature(
                    SimpleLocation(int(comp_start) - 1, int(comp_end), int(comp_strand)),
                    type='CDS',
                    qualifiers={'gene':[comp_gene]},
                )
                # missing ltag means it was not originally annotatable as a standalone gene,
                # so mark it pseudo
                if comp_ltag == '.':
                    comp.qualifiers['pseudo'] = ['']
                else:
                    comp.qualifiers['locus_tag'] = [comp_ltag]
                fusions_by_ltag[ltag]['components'].append(comp)
    return fusions_by_ltag
