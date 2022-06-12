# Release Notes

## Development

### Bugs fixed
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
