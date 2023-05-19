from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

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
            source_alignment=None,
            # pseudo-related attributes
            d3=None, # divisible by 3
            rcs=None, # has reference-corresponding start
            rce=None, # has reference-corresponding end
            vs=None,  # begins with a valid start codon
            ve=None, # ends with a valid stop codon
            bok=None, # BLAST hit ok (meeting thresholds)
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
        self.source_alignment = None

        self.d3 = d3
        self.rcs = rcs
        self.rce = rce
        self.vs = vs
        self.ve = ve
        self.bok = bok
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
