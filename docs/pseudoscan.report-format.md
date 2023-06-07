- locus_tag:
Locus tag of the annotated feature.
In the hybran pipeline, pseudoscan is run before final locus tags have been determined, so the locus tags in those reports will still correspond to those assigned by RATT and the ab initio annotation.
- gene_name:
Assigned gene name (lifted over from reference annotation)
- div_by_3:
Whether the feature sequence is divisible by three.
- valid_start_codon:
Whether the feature sequence begins with a valid start codon according to the detected genetic code table.
- valid_stop_codon:
Whether the feature sequence ends with a valid stop codon according to the detected genetic code table and has no internal stops.
- ref_corr_start:
Whether the feature sequence's beginning aligns to the reference sequence's beginning.
Start refers to the part of the sequence containing the start codon, even for genes on the minus strand.
- ref_corr_end:
Whether the feature sequence's end aligns to the reference sequence's end.
End here refers to the part of the sequence containing the stop codon, even for genes on the minus strand.
- blast_ok:
Whether the feature sequence has a passing BLASTP hit based on the thresholds configured with `--seq-ident` and `--seq-covg`.
- pseudo:
Whether the feature has been called `pseudo`.
- evidence_codes:
Summary of the reason(s) for the `pseudo` determination.
If multiple evidence codes apply to a single feature, they will be comma-delimited.
See below for a description of the possible evidence codes.

Apart from `locus_tag`, `gene_name`, and `evidence_codes`, the column values are `0` (false), `1` (true), or `.` (not determined).

# Evidence Codes
