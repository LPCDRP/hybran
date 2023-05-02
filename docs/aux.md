# Auxiliary Programs

Hybran comes with associated functionality that can be run independently:

`hybran standardize`
: Remove generic gene names.


# `hybran standardize`

`hybran` leaves no gene without a name and only assigns names from genes in your provided references.
Genes that only had *ab initio*-predicted gene names have them relegated to the `gene_synonym` field to avoid propagation during the hybran pipeline.
The result is that if no reference gene name could be assigned, or if a gene is duplicated in the reference genome, a generic name is generated and assigned.

Generic gene names allow easy matching of hybran-identified homologs within and between your samples.
For the purposes of uploading your annotation to a public database, however, these names should not be kept.
`hybran standardize` will remove all these generic names and replace them, wherever possible, with a reference name (for reference duplicates that had been grouped under a generic name) or the *ab initio* predicted name.
