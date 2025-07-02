from copy import copy, deepcopy
import functools

from Bio import (
    GenBank,
    SeqIO,
)
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import translate as super_translate
from Bio.SeqIO import InsdcIO
from Bio.SeqFeature import (
    SeqFeature,
    SimpleLocation,
    CompoundLocation,
    LocationParserError,
)


loc_props = [
    'alignment',
    'd3',
    'rcs',
    'rce',
    'vs',
    've',
    'bok',
    'alts',
    'alte',
    'de',
    'ps_evid',
    'transl_except',
]

class FeatureProperties():
    def __init__(
            self,
            # coordcheck-related attributes
            location=None,
            alignment=None,
            # pseudo-related attributes
            d3=None, # divisible by 3
            rcs=None, # has reference-corresponding start
            rce=None, # has reference-corresponding end
            vs=None,  # begins with a valid start codon
            ve=None, # ends with a valid stop codon
            bok=None, # BLAST hit ok (meeting thresholds)
            alts=None, # has alternative start
            alte=None, # has alternative stop
            de=None, # has delayed end
            ps_evid=None, # list of pseudoscan evidence codes
            transl_except=None, # string of valid 'transl_except' qualifier
    ):
        self.location = location
        self.alignment = alignment

        self.d3 = d3
        self.rcs = rcs
        self.rce = rce
        self.vs = vs
        self.ve = ve
        self.bok = bok
        self.alts = alts
        self.alte = alte
        self.de = de
        if transl_except is None:
            self.transl_except = []
        else:
            self.transl_except = transl_except

        if ps_evid is None:
            self.ps_evid = []
        else:
            self.ps_evid = ps_evid

    def __repr__(self):
        return f"""
================================================================================
Location: {self.location}


A L I G N M E N T
================================================================================
{self.alignment}

================================================================================
Divisible by 3:           {self.d3}
Ref-Corresponding Start:  {self.rcs}
Ref-Corresponding End:    {self.rce}
Valid Start Codon:        {self.vs}
Valid Stop Codon:         {self.ve}
BLAST hit OK:             {self.bok}
Alternative Start Codon:  {self.alts}
Alternative Stop Codon:   {self.alte}
Delayed End:              {self.de}
--------------------------------------------------------------------------------
Pseudoscan Evidence Code: {','.join(self.ps_evid)}
================================================================================
"""


class AutarkicSeqFeature(SeqFeature):
    def __init__(
            self,
            #
            # original SeqFeature attributes
            #
            location=None,
            type='',
            id='<unknown_id>',
            qualifiers=None,
            #
            # our custom attributes
            #
            references=None, # dictionary of `ref` id strings to Seq/SeqRecord objects
            prev_feature=None,
            next_feature=None,
            # similar to prev/next, but implying prev/next on the same strand
            upstream_feature=None,
            downstream_feature=None,
            ###
            source=None, # reference annotation from which this is lifted-over, if applicable
            og=None,
            corr=None,
            corr_possible=None, # whether a coordinate correction was found
            corr_accepted=None, # whether a correction was found not rejected
            ###
            fusion_type=None, # "whole", "partial", or None (if not a fusion at all)
            fusion_components=None, # list of SeqFeatures if fusion_type not None
    ):
        super().__init__(
            location=location,
            type=type,
            id=id,
            qualifiers=qualifiers,
        )

        self.references = references

        self.prev_feature = prev_feature
        self.next_feature = next_feature
        self.upstream_feature = upstream_feature
        self.downstream_feature = downstream_feature

        self.source = source

        if og is None:
            self.og = FeatureProperties()
        else:
            self.og = og
        if corr is None:
            self.corr = FeatureProperties()
        else:
            self.corr = corr
        self.corr_possible = corr_possible
        self.corr_accepted = corr_accepted

        self.fusion_type = fusion_type
        if fusion_components and self.fusion_type is not None:
            self.fusion_components = fusion_components
        else:
            self.fusion_components = []


    # https://stackoverflow.com/a/141777
    @classmethod
    def fromSeqFeature(cls, feature, *args, **kwargs):
        return cls(*args, **vars(feature), **kwargs)

    def extract(self, parent_sequence=None, references=None):
        """
        Make extract() default to using the included parent sequence.
        """
        if references is None and self.references is not None:
            references=self.references
        return super().extract(parent_sequence, references)

    def __getattr__(self, name):
        if name in loc_props:
            if self.corr_accepted:
                props = object.__getattribute__(self, f'corr')
            else:
                props = object.__getattribute__(self, f'og')
            return object.__getattribute__(props, name)
        else:
           return  object.__getattribute__(self, name)

    def __setattr__(self, name, val):
        if name in loc_props:
            if self.corr_accepted:
                props = object.__getattribute__(self, f'corr')
            else:
                props = object.__getattribute__(self, f'og')
            object.__setattr__(props, name, val)
        else:
            object.__setattr__(self, name, val)

    def __eq__(self, other):
        if (
                self.fusion_components and other.fusion_components
                and len(self.fusion_components) == len(other.fusion_components)
        ):
            for i in range(len(self.fusion_components)):
                a = self.fusion_components[i]
                b = other.fusion_components[i]

                a_interesting_qualifiers = {k:v for k, v in a.qualifiers.items() if k in ['locus_tag', 'gene']}
                b_interesting_qualifiers = {k:v for k, v in b.qualifiers.items() if k in ['locus_tag', 'gene']}
                if (
                        a_interesting_qualifiers != b_interesting_qualifiers
                        or a.location != b.location
                ):
                    return False
        elif self.fusion_components != other.fusion_components:
            return False
        return (
            super(AutarkicSeqFeature, self).__eq__(other)
            and self.fusion_type == other.fusion_type
        )

    def __repr__(self):
        fusion_info = ''
        if self.fusion_type:
            fusion_info = [
                f"\nFusion type: {self.fusion_type}",
                "Components:"
            ]
            for c in self.fusion_components:
                component_info = []
                if 'locus_tag' in c.qualifiers:
                    component_info.append(c.qualifiers['locus_tag'][0])
                if 'gene' in c.qualifiers:
                    component_info.append(c.qualifiers['gene'][0])
                component_info.append(str(c.location))
                fusion_info.append(f"{'|'.join(component_info)}")
        return f"""
{self.location}{'\n'.join(fusion_info)}
Qualifiers:
{self.qualifiers}
"""

def translate(
        sequence,
        table='Standard',
        stop_symbol='*',
        to_stop=False,
        cds=True,
        gap=None,
):
    """
    Override the default behavior of Bio.Seq.translate (imported here as super_translate) to prioritize using cds=True unless we can't (due to internal stops and the exception that Biopython throws as a result).
    Ideally, and as a future TODO, we'd like to combine the functionality of `to_stop` and `cds` so that we can have handling of alternative start codons regardless of internal stops.

    All parameters and returns are the same as Bio.Seq.translate.
    The difference is that cds defaults to True and does not conflict with to_stop.
    """

    st = functools.partial(
        super_translate,
        sequence,
        table=table,
        stop_symbol=stop_symbol,
        gap=gap,
    )

    if cds:
        try:
            tr = st(cds=True)
        except TranslationError:
            # TODO: temporarily reset the location of the sequence to exclude everything past the first stop, then call st(cds=True) again.
            tr = st(to_stop=to_stop)
    else:
        tr = st(cds=cds, to_stop=to_stop)

    return tr

# keep track of the original Location.fromstring function before we mask it
super_fromstring = GenBank.Location.fromstring

def locfromstring(text, length=None, circular=False, stranded=True):
    try:
        return super_fromstring(text, length, circular, stranded)
    except:
        # This error will lead Biopython to set the location attribute to None.
        # We'll check for that later and remove it during validation.
        raise LocationParserError from None

class HybGenbankIterator(InsdcIO.GenBankIterator):
    def parse(self, handle):
        records = HybGenBankScanner(debug=0).parse_records(handle)
        return records

class HybGenBankScanner(GenBank.Scanner.GenBankScanner):
    # mask the function used by consumer.location() (_FeatureConsumer)
    GenBank.Location.fromstring = locfromstring

    # - override `parse_features()` to remove None's from the final list of features
    def parse_features(self, skip=False):
        features = super().parse_features(skip=skip)
        return [f for f in features if f is not None]

    # - override `parse_feature()` to return None when encountering any exception
    def parse_feature(self, feature_key, lines):
        try:
            return super().parse_feature(feature_key, lines)
        except:
            return None

class HybGenBankWriter(InsdcIO.GenBankWriter):
    # https://github.com/biopython/biopython/blob/master/Bio/SeqIO/InsdcIO.py
    def _write_feature(self, feature, record_length):
        """
        intercept the .ref attributes to prevent locations from being printed
        with f"{ref}:" sort of prefixed to it.
        """
        # our AutarkicSeqFeature stores the contig associated with the feature, so
        # regular copy here instead of deepcopy to avoid spiking memory usage.
        temp_feature = copy(feature)
        # But we need a deepcopy of the location because CompoundLocations have
        # component SimpleLocations and we don't want to change the originals
        temp_loc = deepcopy(feature.location)
        temp_feature.location = temp_loc
        for part in temp_feature.location.parts:
            part.ref = None
        return super()._write_feature(temp_feature, record_length)

SeqIO._FormatToWriter['genbank'] = SeqIO._FormatToWriter['gb'] = HybGenBankWriter
SeqIO._FormatToIterator['genbank'] = SeqIO._FormatToIterator['gb'] = HybGenbankIterator
