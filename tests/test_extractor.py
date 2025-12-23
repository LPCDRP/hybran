from Bio.SeqFeature import (
    CompoundLocation,
    SimpleLocation,
    ExactPosition,
)

from hybran import extractor
from hybran.util import keydefaultdict

from .data_features import *


def test_ref_fuse():
    ref_annotation = keydefaultdict(extractor.ref_fuse)
    ref_annotation.update(ref_features['H37Rv'])
    ref_id = 'H37Rv.NC_000962.3'

    assert ref_annotation[f'{ref_id}::{ref_id}@@@PE_PGRS50::PE_PGRS49'].location == CompoundLocation([
        SimpleLocation(ExactPosition(3738157), ExactPosition(3742774), strand=-1, ref=ref_id),
        SimpleLocation(ExactPosition(3736983), ExactPosition(3738000), strand=-1, ref=ref_id)
    ], 'join')
