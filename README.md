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
│   │   ├── hybran_coord_corrections.tsv
│   │   ├── merged_genes.gbk
│   │   ├── prokka_unused.tsv
│   │   └── ratt_unused.tsv
|   ├── ratt/
│   │   └── ...
|   ├── prokka/
│   │   ├── hybran_coord_corrections.tsv
│   │   └── ...
├── <i>sampleN</i>/
│   └── ...
│
├──deduped-refs/ (if run with <i>--dedupe-references</i>)
│   ├── duplicates.tsv
│   ├── <i>reference1</i>.gbk
│   ├── ...
│   └── <i>referenceN</i>.gbk
├── prodigal-test/
│   ├── incorrect_starts.txt
│   └── reference_prodigal
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

##### `deduped-refs/duplicates.tsv`

If `hybran` was called with `--dedupe-references`, revised reference annotations will be generated in the `deduped-refs` directory.
These annotations differ from the original in that each set of duplicate genes is assigned a single name used for all instances.
The original name is retained as a `gene_synonym` qualifier in the annotation file.
The file `duplicates.tsv` will list duplicate genes found in the reference annotations and the name they were assigned.
Columns in this file are

* reference name
* reference locus tag
* reference gene name
* unified name

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

##### `annomerge/*/{ratt,prokka}_unused.tsv`

Two columns: the locus tag and the reason why it was excluded from the final annotation.


##### `*/*/hybran_coord_corrections.tsv`


- reference_locus_tag:
The gene from your reference annotation corresponding to the annotated feature.
- locus_tag:
Annotated feature from the source indicated by the parent directory.
- original_start
- original_end
- original_strand
- start
- end
- strand
- start_shift
The difference between `start` and `original_start`.
End coordinates will not have been adjusted.

##### `prodigal-test/incorrect_starts.txt`

The results of the comparison of the reference annotation(s) to Prodigal's annotation of them.
Genes with differing start positions are looked for in each sample and their coordinates checked more closely against the reference annotation.

- reference_locus_tag
- length_according_to_prodigal
- start_change


## Citation
[Manuscript submitted]
