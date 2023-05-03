# Commands and Auxiliary Programs

## Core Subcommands
The core elements of the Hybran pipeline can be run separately, but they are being added incrementally:

`hybran onegene`
: Unify names of gene duplicates.

### `hybran onegene`

`onegene` is used to process the reference annotations to avoid inconsistent matching to different instances of (near-)identical genes.
It solves the problem by assigning a single generic name to all instances, while storing the original name in each feature's `gene_synonym` qualifier.
The unified generic name propagates to the input samples' final annotations, but if this is not desired, you can use `hybran standardize` on the output annotations to restore the original names.


## Auxiliary Programs

Hybran comes with associated functionality that can be run independently.
These are either not part of the pipeline at all or not a component in their own right.
These tools are:

`hybran standardize`
: Remove generic gene names.


### `hybran standardize`

`hybran` leaves no gene without a name and only assigns names from genes in your provided references.
Genes that only had *ab initio*-predicted gene names have them relegated to the `gene_synonym` field to avoid propagation during the hybran pipeline.
The result is that if no reference gene name could be assigned, or if a gene is duplicated in the reference genome, a generic name is generated and assigned.

Generic gene names allow easy matching of hybran-identified homologs within and between your samples.
For the purposes of uploading your annotation to a public database, however, these names should not be kept.
`hybran standardize` will remove all these generic names and replace them, wherever possible, with a reference name (for reference duplicates that had been grouped under a generic name) or the *ab initio* predicted name.
