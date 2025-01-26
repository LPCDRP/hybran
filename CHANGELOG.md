# Release Notes

## Development

### Enhancements
* Produce gff files in addition to gbk for onegene, standardize, and postprocessing
* Write a log file to the output folder (#73)

### Bugs fixed
* Fixed verbose and quiet run modes.
* The criteria for accepting coordinate corrections now also requires valid start codons. (#71)
* Fixed off by 1 error when reporting the location of translational exceptions (selenocysteine proteins). (#72)
* Lifted memory restrictions on cd-hit (#74)
* Fixed crash that occurs when only one gene component of a fusion doesn't have a product field. (#75)

## [Version 1.8](https://gitlab.com/LPCDRP/hybran/-/tags/1.8)

* hybran now depends on `networkx` and `intervaltree`

### Enhancements
* New subcommand:
  - `hybran compare`: compare two annotations of the same genome
* Rewritten comparison engine based on intervaltrees.
  This fixes some edge cases, such as pairs of in-frame overlapping genes with
  an unrelated gene between them, as well as unannotated regions at contig edges.
  The code is modular, available to run in the above-mentioned new command, and used in annomerge.
* Enhancements to `hybran standardize`:
  - New option `-r`/`--ref-names-only`.
  - Simplified usage: If passing in a hybran output directory, no other arguments will be necessary.
* Annotations now include pseudoscan evidence codes as note fields.
These were previously only found in the pseudoscan report files.

### Bugs fixed
* Annomerge no longer produces redundant CDS entries (#69)
* Synonymous SNPs in start codons do not automatically invalidate reference-correspondence.
* pseudoscan: Fixed typo potentially affecting calls where the reference is pseudo (#68).
* fusionfisher:
  - Identify another signature of hybrid fusion genes (previously identified as misannotations)
  - Update FeatureProperties for fusion genes based on their components. Fusions no longer categorically pseudo.
  - Eliminated more false positive fusion calls by accounting for another misannotation signature.
* Coordinate corrections no longer collapse compound locations ("joins") into single intervals (#70).
* Include handling of `transl_except` qualifiers for amino acids introduced through translational recoding.
* Restored compatibility with biopython >=1.82 and applied more preventitive maintenance in response to deprecation warnings.
* Translations, in particular for RATT-transferred annotations, fixed to account for alternative methionine start codons.
* Recognize .gbff as a genbank format file extension.
* Write out all Prokka postprocessed contigs' annotations instead of just the last one.

## [Version 1.7.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.7.1)

### Bugs fixed
* Fixed issue with `onegene` where unnamed reference genes did not get unified to use the locus tag of the primary reference as a gene name, instead unnecessarily assigning a generic name.

## [Version 1.7](https://gitlab.com/LPCDRP/hybran/-/tags/1.7)

### Enhancements
* New subcommands:
  - `hybran standardize`: remove generic gene names in the final annotations.
  - `hybran onegene`: unify names of highly conserved gene copies
* New option `-s`/`--organism` to set organism name in the final genbank header (#58).
* Intercept malformed features in RATT output to allow continuity of the pipeline (#64).
*  The `--dedupe-references` option is now deprecated, being made a core part of the pipeline with `onegene` since the generic names that it introduces to collapse reference paralogs can be undone by `hybran standardize` if desired afterwards.
* Reference unification now uses a different generic prefix than that used for unnamed genes to allow differentiating between them.
`hybran standardize` handles these as well.
* Full support for multiple reference annotations.
* Postprocessed versions of RATT and ab initio annotations, with associated reports, are now saved.
* Substantial speedups due to parallelization of all postprocessing logic.
* Logging and reporting:
  - Comprehensive reporting of pseudoscan results.
  - Reorganized the invalid and rejected features reports.
  - Removed spaces from column names of novelty report.

### Bugs fixed
* reference unification (formerly activated by the `--dedupe-references` option):
  - fixed for multiple references (#55)
  - enabled checkpointing on this step in the full pipeline.
* pseudoscan: now identifies delayed stop codons when start coordinate correction changes the reading frame without introducing an internal stop.
* coord_check:
  - Fixed issue in coord_check's reporting status when not attempting correction.
  - Coordinate correction sometimes resulted in genes with no stop codons.
    coord_check now checks for this and extends the ORF to the next in-frame stop codon to make a proper correction.
  These address the root cause of the problem for which the temporary measure from v1.6.1 was taken.
  That temporary measure has been removed.
* Improved detection of gene fusions when one of the components is derived from a reference pseudogene.
* fissionfuser: fixed issue where apparent complementary fragments are combined despite one of the copies being non-pseudo (#66).
Although the issue severity was mitigated in that the combined annotation would be rejected in favor of RATT's, the second gene would have been lost.
* annomerge: reject the `source` feature from RATT.
Some final annotations contained two of them: one from RATT and one from Prokka.
* Fixed issues that occur when sequence IDs contain "|" character (#62).
* Fixed handling of situations where either RATT or Prokka find no annotations.
* Fixed problem with redundant fusion gene name components when detected using both RATT and Prokka.
* Improved detection of gene fusions due to adjustments of alignment internal gap extension penalty and refined delayed-stop calling criteria.

### Housekeeping
* Changed default generic ORF prefix to "HYBRA" for greater clarity.
* set default output directory to current directory.
* Now ignoring warnings when generating translations for pseudogenes that aren't multiples of three.
* Enabled setting any of RATT's configured transfer types and fixed names for *.global parameter sets.

## [Version 1.6.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.6.1)

### Bugs fixed
* Fixed scenario where hybran crashes when coordinate correction matches to an adjacent locus and attempts correction (#61).
* Added a temporary measure to inclusion criteria to penalize CDSs lacking stop codons (inadequately postprocessed) (#60).

## [Version 1.6](https://gitlab.com/LPCDRP/hybran/-/tags/1.6)

### Enhancements
* Massive streamlining of the pipeline. Reworked components into new subsystems:
  - `pseudoscan`: identification of anomalous copies of reference genes using new criteria independent of alignment coverage. (#50, #56, and #59)
  - `fissionfuser` (formerly `process_split_genes()`): improved detection and combining of gene fragments (that ab initio annotations tend to produce) into a single record.
  - `fusionfisher`: detection of gene fusion events and putative misannotations.
  - `thunderdome`: more aggressive conflict resolution between RATT and ab initio annotations.
* Output GFF files no longer include the genome sequence.

### Bugs fixed
* Fixed handling of conflicting annotations that are differently named (#57).
* Reimplemented coordinate correction and applied to ab initio ORFs as part of `pseudoscan`.
  This resolves many instances of false `pseudo` CDSs ab initio that were due simply to
  incorrect start coordinate predictions spuriously shortening the genes.
* Fixed handling of compound intervals in reference annotations (#46, #47)
* Resolved issues involving reference annotations with multiple contigs/chromosomes (#48)
* Fixed issue with some gene name assignments being dropped later in the pipeline due to some obsolete code (#43).
* More comprehensive tracking of RATT/ab initio overlaps and conflicts (#49).
* Checking in-frame overlaps with pseudo ORFs containing internal stop codons
* Revamped postprocessing of RATT-introduced compound intervals (#44, #45)
* Updated inclusion criteria for special handling of pseudo ORFs (#42)

## [Version 1.5.2](https://gitlab.com/LPCDRP/hybran/-/tags/1.5.2)

### Bugs fixed
* Made a consistent non-CDS policy for RATT: Take everything except rRNA and tRNA (#22)
* Clarified some rejection reasons for RATT/ab initio features.
* Fixed representation of blast results for CDSs when there aren't any hits
* Fixed issue with RATT handoff if sample/contig names contain `.` or `|`.
* Fixed issue with logging merged genes.

## [Version 1.5.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.5.1)

### Bugs fixed
* Fixed newly introduced issue with rejecting RATT annotations.

### Housekeeping
* Removed some unused code.

## [Version 1.5](https://gitlab.com/LPCDRP/hybran/-/tags/1.5)

### Bugs fixed
* Prevented overlapping RATT-transferred annotations from automatically being handled as conflicts, leading one of the two to be discarded (#39)
* Corrected distinguishing of ab initio vs reference-transferred annotation in final conflict resolution step (#38)
* Added logging of some missed cases of annotation rejections (for {ratt,prokka}_unused.tsv)
* Fixed the reference gene <=> locus_tag mapping dictionaries used in annomerge (#41)
* Made sure to track and process all ab initio annotations that overlap RATT-transferred CDSs (#40).
* Fixed handling of multi-fasta inputs (#35)

### Enhancements
* Streamlined Prokka workflow (#33)
* Parallelized BLASTing to reference genes in annomerge.
* Set sequence names in the output annotation files.

### Housekeeping
* Added exit status checks so pipeline fails as early as possible when things go wrong.
* Switched default evalue to Prokka's current setting of 1e-9 (was 1e-6)
* Slightly streamlined the RATT / Prokka comparison workflow
* Added more unit tests

## [Version 1.4.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.4.1)

### Bugs fixed
* Genes split into multiple adjacent fragments used to have a single /gene record but multiple CDS records with the same locus tag.
  For INSDC compliance, they now only have a single CDS record as well.
* Removed /translation fields for /pseudo CDSs.

## [Version 1.4](https://gitlab.com/LPCDRP/hybran/-/tags/1.4)

### Enhancements
* Generalized for any prokaryote.
    - Genetic code and taxonomy ID detected from reference annotation.
      RATT configuration is now automatically generated based on the
      detected genetic code, so a configuration file is no longer bundled.
    - Now using "ORF" prefix for generic genes rather than "MTB" by default.
      Option `--orf-prefix` added for customizability.
* Removed checking for dnaA as the first gene at the first base position.
* Made eggNOG-mapper step optional.
* Gene fragments are now identified using the corresponding reference gene names, but are distinguished with a /pseudo tag.
* RATT and (some) Prokka options are now under user control.

### Bugs fixed
* Account for translationless CDSs that are labeled with the 'pseudo' qualifier instead of 'pseudogene'
* Allow input fasta files with alternative standard extensions.
* Fixed handling of reference annotations that may not have /gene qualifiers for all annotations.
* Fixed handling of input genome when it's the same as the reference.
* Set proper field from which to draw eggNOG-mapper annotations.
* Uniform locus tags are now assigned for every sample.
* Better identification of reference and unnamed genes in processing of clusters.

## [Version 1.3.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.3.1)

### Enhancements
* Hybran version now recorded in the genbank annotation header.

### Bugs fixed
* Updated the Prokka reference proteome generation format to enable Prokka to set gene names and product fields rather than leaving it to the final clustering step.
* Fixed installation location of resource file.

## [Version 1.3.0](https://gitlab.com/LPCDRP/hybran/-/tags/1.3.0)

### Enhancements
* Added `--dedupe-references` option to assign a single generic gene name to duplicate genes in the provided reference annotations.
* Sequence identity and alignment coverage thresholds are no longer applied to RATT-transferred annotations by default. (#28)
  The original behavior can be restored by passing the new `--filter-ratt` option.
* reference annotations can now be passed as individual file names or file of file names, in addition to a directory name (#21)

### Bugs fixed
* eggnog-mapper step no longer gets skipped (#23)
* alignment query coverage threshold is now applied directly in Prokka (#24)
* Dropped criterion of excluding hypothetical genes from Prokka-no-reference (#30)
* Fixed calculation of query and reference alignment coverage (#27).
* Corrected selection of top blastp hits for the one-to-one and one-to-many searches.
  When there were multiple hits in these cases, only the last one output by BLAST was being retained, which actually corresponds to the worst hit (by e-value).
  We now retain only the first hit.


## [Version 1.2.0](https://gitlab.com/LPCDRP/hybran/-/tags/1.2.0)

### Enhancements
* Add thorough logging of gene annotations merged, rejected, and newly-named. (#19)
* Better tolerance of directory name inputs (#18)
* New option `-c`/`--coverage-threshold` for tuning gene matching.
  The sequence identity threshold is now taken through `-i`/`--identity-threshold`.
* Incorporate a RATT configuration file for using codons from translation table 11 (#17)

### Bugs fixed
- use provided identity/coverage thresholds for all instances of
BLAST, CD-HIT (#10)
- ensure that the final gff output file gets updated (#15)
- (#14)


## [Version 1.1.1](https://gitlab.com/LPCDRP/hybran/-/tags/1.1.1)
* Fixed issue with writing of merged_genes.gbk that sometimes
  caused hybran to crash during a run.


## [Version 1.1.0](https://gitlab.com/LPCDRP/hybran/-/tags/1.1.0)

* Renamed to Hybran
* Migrated to Python 3
* Removed limitation of 30 references
* Allowed sequence identity threshold to be user-defined
* Proper handling of temporary files
* Fixed issue preventing clustering step from running

## [Version 1.0](https://gitlab.com/LPCDRP/hybran/-/releases#v1.0)

Initial release of AnnoTUB
