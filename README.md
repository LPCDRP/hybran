# AnnoTUB (Annotate TUBerculosis)

Annotation pipeline that provides genomic annotation 
for *Mycobacterium tuberculosis* de novo assembled genomes

### Dependencies
* [RATT](http://ratt.sourceforge.net/)
* [Prokka](https://github.com/tseemann/prokka)
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [python2.7](https://www.python.org/downloads/release/python-2715/)
* [Biopython](https://biopython.org/wiki/Download)
* [numpy](https://sourceforge.net/projects/numpy/files/NumPy/)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [MCL](https://github.com/JohannesBuchner/mcl)
* [NCBI-BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)

### How to run
This can be executed on one or many genomes. The more reference
annotations included, the more accurate the annotation will be 
and less ambiguity will exist for the target genomes.
```
./annotub                                        \
    --genomes /dir/to/FASTAs                   \
    --references /dir/to/reference/annotations \
    --emapper /dir/to/eggnog-mapper            \
    --output ./                                \
    --nproc 2
```

### How to cite
Manuscript submitted
