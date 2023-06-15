import Bio
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord

from .bio import AutarkicSeqFeature, SeqIO
from . import designator
from .extractor import get_ltag, get_gene

def log_feature_fate(feature, logfile, superior=None, evid=".", remark="."):
    """
    General-purpose logging function to print out a gene's information and a comment
    :param feature: A SeqFeature object
    :param logfile: An open filehandle
    :param superior: SeqFeature representing the superior annotation displacing feature.
    :param evid: (str) evidence code
    """
    if superior is not None:
        superior_locus_tag = get_ltag(superior)
        superior_gene_name = get_gene(superior)
    else:
        superior_locus_tag = superior_gene_name = '.'

    locus_tag = get_ltag(feature)
    gene_name = get_gene(feature)

    print('\t'.join([
        locus_tag,
        gene_name,
        evid,
        remark,
        superior_locus_tag,
        superior_gene_name,
    ]), file=logfile)

def log_coord_correction(feature, logfile):
    """
    This function is used to log gene information when coord_check() fixes a start/stop postition for individual features.
    It is called by log_coord_corrections().
    :param feature: AutarkicSeqFeature object
    :param logfile: An open filehandle
    """
    locus_tag = feature.qualifiers['locus_tag'][0]
    gene_name = feature.qualifiers['gene'][0]
    strand = str(feature.strand)
    og_start = (int(feature.og.location.start) + 1)
    og_end = (int(feature.og.location.end))
    new_start = (int(feature.corr.location.start) + 1)
    new_end = (int(feature.corr.location.end))
    start_fixed = str(og_start != new_start).lower()
    stop_fixed = str(og_end != new_end).lower()
    gene_length_ratio = f"{(og_end - og_start)/(new_end - new_start):.3f}"

    if feature.strand == -1:
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
        'div_by_3',
        'valid_start_codon',
        'valid_stop_codon',
        'ref_corr_start',
        'ref_corr_end',
        'blast_ok',
        'pseudo',
        'evidence_codes',
    ]
    # represent boolean values as ints but account for N/A (None)
    nacast = lambda _: int(_) if _ is not None else "."
    print('\t'.join(header), file=logfile)
    for contig in features_by_contig_dict:
        for feature in features_by_contig_dict[contig]:
            if feature.ps_evid:
                line = [
                    feature.qualifiers['locus_tag'][0],
                    feature.qualifiers['gene'][0],
                    nacast(feature.d3),
                    nacast(feature.vs),
                    nacast(feature.ve),
                    nacast(feature.rcs),
                    nacast(feature.rce),
                    nacast(feature.bok),
                    nacast(designator.is_pseudo(feature.qualifiers)),
                    ','.join(feature.ps_evid),
                ]
                print('\t'.join(str(v) for v in line), file=logfile)
