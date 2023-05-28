from copy import copy

from Bio import SeqIO
from Bio.SeqIO import InsdcIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation


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
    ):
        self.location = location
        self.alignment = alignment

        self.d3 = d3
        self.rcs = rcs
        self.rce = rce
        self.vs = vs
        self.ve = ve
        self.bok = bok


class AutarkicSeqFeature(SeqFeature):
    def __init__(
            self,
            #
            # original SeqFeature attributes
            #
            location=None,
            type='',
            location_operator='',
            strand=None,
            id='<unknown_id>',
            qualifiers=None,
            sub_features=None,
            ref=None,
            ref_db=None,
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
            og=FeatureProperties(),
            corr=FeatureProperties(),
            corr_possible=None, # whether a coordinate correction was found
            corr_accepted=None, # whether a correction was found not rejected
            ###
            is_fusion=None, # is a fusion of two or more genes
    ):
        super().__init__(
            location=location,
            type=type,
            location_operator=location_operator,
            strand=strand,
            id=id,
            qualifiers=qualifiers,
            sub_features=sub_features,
            ref=ref,
            ref_db=ref_db,
        )

        self.references = references

        self.prev_feature = prev_feature
        self.next_feature = next_feature
        self.upstream_feature = upstream_feature
        self.downstream_feature = downstream_feature

        self.source = source

        self.og = og
        self.corr = corr
        self.corr_possible = corr_possible
        self.corr_accepted = corr_accepted

        self.is_fusion = is_fusion


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



class HybGenBankWriter(InsdcIO.GenBankWriter):
    # https://github.com/biopython/biopython/blob/master/Bio/SeqIO/InsdcIO.py
    def _write_feature(self, feature, record_length):
        """
        intercept the .ref attributes to prevent locations from being printed
        with f"{ref}:" sort of prefixed to it.
        """
        temp_feature = copy(feature)
        temp_loc = copy(feature.location)
        temp_feature.location = temp_loc
        for part in temp_feature.location.parts:
            part.ref = None
        return super()._write_feature(temp_feature, record_length)

SeqIO._FormatToWriter['genbank'] = SeqIO._FormatToWriter['gb'] = HybGenBankWriter
