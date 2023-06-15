# Hybran

Hybran is a hybrid reference-based and ab initio genome annotation pipeline for prokaryotic genomes.
It uses the [Rapid Annotation Transfer Tool (RATT)](http://ratt.sourceforge.net) to transfer as many annotations as possible from your reference genome annotation based on conserved synteny between the nucleotide genome sequences.
Hybran then supplements unannotated regions with *ab initio* predictions from [Prokka](https://github.com/tseemann/prokka).
Then, all coding sequence annotations are clustered together and additional reference gene names are assigned based on amino acid sequence identity and alignment coverage.

## Execution

This can be executed on one or many genomes. The more reference
annotations included, the more accurate the annotation will be 
and less ambiguity will exist for the target genomes. Input can
be a FASTA, a list of FASTAs (space-separated), a directory containing
FASTAs, or a File Of FileNames (FOFN) of FASTAs.
```
hybran                                                                          \
    --genomes /dir/to/FASTAs | in.fasta [in2.fasta in3.fasta ...] | fastas.fofn  \
    --references /dir/to/reference/annotation(s)                                 \
    --output ./                                                                  \
    --nproc 2
```

## Output

Final annotations are created in Genbank and GFF formats in the output directory.
The output directory also contains intermediate files and informative logs:

<pre><code>
outdir/
├── <i>sample1</i>/
│   ├── annomerge/
│   │   ├── <i>sample1</i>.gbk
│   │   ├── <i>sample1</i>.gff
│   │   ├── coord_corrections.tsv
│   │   ├── prokka_unused.tsv
│   │   ├── pseudoscan_report.tsv
│   │   └── ratt_unused.tsv
|   ├── ratt/
│   │   └── ...
|   ├── ratt-postprocessed/
│   │   ├── <i>sample1</i>.*.final.gbk
│   │   ├── coord_corrections.tsv
│   │   ├── invalid_features.tsv
│   │   └── pseudoscan_report.tsv
|   ├── prokka/
│   │   └── ...
|   ├── prokka-postprocessed/
│   │   ├── <i>sample1</i>.gbk
│   │   ├── coord_corrections.tsv
│   │   ├── invalid_features.tsv
│   │   └── pseudoscan_report.tsv
├── <i>sampleN</i>/
│   └── ...
│
├──unified-refs/
│   ├── unifications.tsv
|   ├── unique_ref_cdss.faa
│   ├── <i>reference1</i>.gbk
│   ├── ...
│   └── <i>referenceN</i>.gbk
├── clustering/
│   ├── multigene_clusters.txt
│   ├── novelty_report.tsv
│   ├── onlyltag_clusters.txt
│   └── singleton_clusters.txt
|
├── <i>sample1</i>.gbk
├── <i>sample1</i>.gff
├── ...
├── ...
├── <i>sampleN</i>.gbk
└── <i>sampleN</i>.gff
</code></pre>

#### Logs and Reports

##### `unified-refs/unifications.tsv`

`hybran` generates revised reference annotations in the `unified-refs` directory.
These annotations differ from the original in that each set of conserved (>=99% amino acid identity and alignment coverage) or duplicated genes is assigned a single name used for all instances.
The original name is retained as a `gene_synonym` qualifier in the annotation file.
The file `unifications.tsv` will list duplicate genes found in the reference annotations and the name they were assigned.
Columns in this file are

* reference name
* reference locus tag
* reference gene name
* unified name

##### `unified-refs/unique_ref_cdss.faa`

A multi-fasta file of the representative amino acid sequences for each unique reference CDS.

##### `clustering/novelty_report.tsv`

Depending on the sequence identity and alignment coverage thresholds used, Hybran will name candidate novel genes.
This novelty report allows you to examine whether these genes are truly unique based on how close they came to meeting the thresholds.

* cluster_type
* candidate_novel_gene
* nearest_ref_match
The top hit among the reference or other candidate novel genes.
* metric
The `nearest_ref_match` is the top hit according to the metric specified in this column.
Its values for all three metrics are shown in the next columns.
* % AA sequence identity
* % subject(ref) coverage
* % query coverage

##### `*/annomerge/{ratt,prokka}_unused.tsv` `*/{ratt,prokka}-postprocessed/invalid_features.tsv`

- locus_tag:
Locus tag of the rejected feature from the source indicated by the file name or parent directory.
- gene_name:
Assigned gene name of the rejected feature (lifted over from reference annotation).
Same as locus tag if none was assigned.
- rival_locus_tag:
Locus tag of the prevailing feature.
- rival_gene_name:
Assigned gene name of the prevailing feature (lifted over from reference annotation).
- evidence_codes:
Summary of the reason for rejecting the feature.
- remark:
A more verbose explanation of the rejection reason.

##### Evidence Codes

###### Evidence Codes Assigned during RATT Postprocessing

- no_coordinates
: RATT sometimes outputs malformed feature locations (see, for example, [RATT#18](https://github.com/ThomasDOtto/ratt/issues/18) and [RATT#19](https://github.com/ThomasDOtto/ratt/issues/19)).
Hybran intercepts these during parsing of the results and sets an empty location to enable continuity of the pipeline.
Since the malformed feature could not be properly parsed, however, there may not be a name to refer to in the log here.
- zero_length
- categorical
: Currently, rRNAs and tRNAs are only taken from the ab initio annotation, so these are categorically rejected from RATT.
- misplaced
- poor_match
: When using `--filter-ratt`, annotations not meeting the blastp thresholds are rejected and have this evidence code applied.

###### Evidence Codes Assigned by fissionfuser

`fissionfuser` is only applied during postprocessing of the ab initio annotations.

- complementary_fragments
- overlapping_inframe
: This scenario arises as a result of postprocessing ab initio annotations.
When a CDS has an internal stop, the ab initio annotation often reports what looks like a tandem duplication.
Start coordinate correction by pseudoscan often extends the downstream fragment to overlap with the upstream fragment and `fissionfuser` identifies this fission event signature.

###### Evidence Codes Assigned by fusionfisher

- redundant_fusion_member
- combined_annotation
- putative_misannotation

###### Evidence Codes Assigned by annomerge

- identical
- identical_non_cds
- shorter
- shorter_pseudogene
- forfeit
: When postprocessed RATT and ab initio annotations are equally valid, RATT is preferred since its name assignment derives from synteny.
- internal_stop
- worse_ref_correspondence
- nonpseudo_vs_pseudo
- hypothetical_vs_real
: When an ab initio annotation for which a name could not be assigned using blastp hits conflicts with a RATT annotation, the ab initio annotation is rejected for this reason.

##### `*/*/coord_corrections.tsv`

- locus_tag:
Locus tag of the annotated feature from the source indicated by the parent directory.
- gene_name:
Assigned gene name (lifted over from reference annotation)
- strand
- og_start
Original start position
- og_end
Original end position
- new_start
Updated start position
- new_end
Updated end position
- fixed_start_codon
Whether the start codon was corrected ('true' or 'false')
- fixed_stop_codon
Whether the stop codon was corrected ('true' or 'false')
- gene_length_diff
The percent difference in gene length between the original and updated locations
- status
Whether the correction was accepted or rejected

For `og_start`, `og_end`, `new_start`, and `new_end`, "start" always corresponds to the low number on the genome and "stop" corresponds to the high number, regardless of strand.
`new_start` and `new_end` are not necessary modified from the original coordinates.
`fixed_start_codon` and `fixed_stop_codon` indicate whether they have changed, but these correspond to the strand-adjusted start and stop positions, hence the reference to codons.


##### `*/*/pseudoscan_report.tsv`

A summary of the characteristics of "interesting" features found by [pseudoscan](pseudoscan.md).
Such features include all genes to which the `pseudo` tag was applied, but also includes non-pseudo genes if they had signatures consistent with a `pseudo` but had a redeeming attribute.

{!pseudoscan.report-format.md!}

## Citation

Elghraoui, A.; Gunasekaran, D.; Ramirez-Busby, S. M.; Bishop, E.; Valafar, F.
Hybran: Hybrid Reference Transfer and Ab Initio Prokaryotic Genome Annotation.
bioRxiv November 10, 2022, p 2022.11.09.515824. doi:[10.1101/2022.11.09.515824](https://doi.org/10.1101/2022.11.09.515824).
