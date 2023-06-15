- locus_tag:
Locus tag of the annotated feature.
In the hybran pipeline, pseudoscan is run before final locus tags have been determined, so the locus tags in those reports will still correspond to those assigned by RATT and the ab initio annotation.
- gene_name:
Assigned gene name (lifted over from reference annotation)
- pseudo:
Whether the feature has been called `pseudo`.
- evidence_codes:
Summary of the reason(s) for the `pseudo` determination.
If multiple evidence codes apply to a single feature, they will be semicolon-delimited.
See below for a description of the possible evidence codes.
- summary:
Semicolon-delimited string of (abbreviated) pseudo attributes and their respective boolean values.
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
Whether the feature sequence has a passing blastp hit based on the thresholds configured with `--seq-ident` and `--seq-covg`.

Apart from `locus_tag`, `gene_name`, and `evidence_codes`, the column values are `0` (false), `1` (true), or `.` (not determined).

# Evidence Codes
- ref_pseudo:
The reference gene is marked `pseudo`.
- alt_start:
The feature has a valid start codon that does NOT correspond with the reference.
- alt_end:
The feature has a valid stop codon that does NOT correspond with the reference.
- noisy_seq:
The feature has a reference-corresponding start and stop, but a poor blastp hit.
- no_rcc:
The feature does not have a reference-corresponding start and stop. The feature also has a poor blastp hit.
- not_div_by_3:
The feature does not have a valid reading frame.
- internal_stop:
The feature contains an internal stop codon.
