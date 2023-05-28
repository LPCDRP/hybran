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
│   │   ├── prokka_unused.tsv
│   │   └── ratt_unused.tsv
|   ├── ratt/
│   │   └── ...
|   ├── ratt-postprocessed/
│   │   ├── <i>sample1</i>.*.final.gbk
│   │   ├── coord_corrections.tsv
│   │   └── invalid_features.tsv
|   ├── prokka/
│   │   └── ...
|   ├── prokka-postprocessed/
│   │   ├── <i>sample1</i>.gbk
│   │   ├── coord_corrections.tsv
│   │   └── invalid_features.tsv
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

Two columns: the locus tag and the reason why it was excluded from the final annotation.


##### `*/{ratt,prokka}-postprocessed/coord_corrections.tsv`

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

For `og_start`, `og_end`, `new_start`, and `new_end`, "start" always corresponds to the low number on the genome and "stop" corresponds to the high number, regardless of strand.
`new_start` and `new_end` are not necessary modified from the original coordinates.
`fixed_start_codon` and `fixed_stop_codon` indicate whether they have changed, but these correspond to the strand-adjusted start and stop positions, hence the reference to codons.

## Citation

Elghraoui, A.; Gunasekaran, D.; Ramirez-Busby, S. M.; Bishop, E.; Valafar, F.
Hybran: Hybrid Reference Transfer and Ab Initio Prokaryotic Genome Annotation.
bioRxiv November 10, 2022, p 2022.11.09.515824. doi:[10.1101/2022.11.09.515824](https://doi.org/10.1101/2022.11.09.515824).
