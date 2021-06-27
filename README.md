# Hybran

Genome annotation pipeline for *Mycobacterium tuberculosis* de novo assembled genomes

### Dependencies
* [RATT](http://ratt.sourceforge.net/)
* [Prokka](https://github.com/tseemann/prokka)
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [Biopython](https://biopython.org/wiki/Download)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [MCL](https://github.com/JohannesBuchner/mcl)
* [NCBI-BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)
    * Required database: bactNOG v4.5.0
    * To download: `download_eggnog_data.py -y bactNOG`
##### All dependencies are available within the `environment.yml` within an anaconda environment

### How to install

```
# Set up the anaconda environment with the dependencies
conda env create -n hybran -f environment.yml
conda activate hybran

# install hybran
python setup.py install
```

### How to run
This can be executed on one or many genomes. The more reference
annotations included, the more accurate the annotation will be 
and less ambiguity will exist for the target genomes. Input can
be a FASTA, a list of FASTAs (space-separated), a directory containing
FASTAs, or a File Of FileNames (FOFN) of FASTAs.
```
hybran                                                                          \
    --genomes /dir/to/FASTAs | in.fasta [in2.fasta in3.fasta ...] | fastas.fofn  \
    --references /dir/to/reference/annotation(s)                                 \
    --eggnog-databases /dir/to/eggnog-mapper-database/data/                      \
    --output ./                                                                  \
    --nproc 2
```

A sample FASTA is available to demo hybran:
```
hybran 
    --genomes sample/
    --references annotations/
    --eggnog-databases eggnog-mapper/data/
    --output ./
    --nproc 4
```
#### Large Number of References

If you are using a large number of references, you may get an error during the annotation transfer step like this:

```
ERROR: mummer and/or mgaps returned non-zero
```

If so, you should recompile mummer using `make CPPFLAGS="-O3 -DSIXTYFOURBITS"` as described at <https://sourceforge.net/p/mummer/mailman/message/34069398/>.

### Output

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
|   ├── prokka-noreference/
│   │   ├── hybran_coord_corrections.tsv
│   │   └── ...
├── <i>sampleN</i>/
│   └── ...
│
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
`prokka_unused.tsv` refers to the Prokka-reference run.


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


### How to cite
[Manuscript submitted]
