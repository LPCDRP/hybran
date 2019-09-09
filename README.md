# AnnoTUB (Annotate TUBerculosis)

Genome annotation pipeline for *Mycobacterium tuberculosis* de novo assembled genomes

### System Requirements
* Linux OS required
* Tested on Debian 8 (jessie)
* Install time ~5 minutes

### Dependencies
* [RATT](http://ratt.sourceforge.net/) >= v18
* [Prokka](https://github.com/tseemann/prokka) >=v1.13
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [python2.7](https://www.python.org/downloads/release/python-2715/)
* [Biopython](https://biopython.org/wiki/Download)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [MCL](https://github.com/JohannesBuchner/mcl)
* [NCBI-BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)

### How to install
```
python2.7 setup.py install --prefix=~/bin/
```

### How to run
This can be executed on one or many genomes. The more reference
annotations included, the more accurate the annotation will be 
and less ambiguity will exist for the target genomes.
```
annotub                                            \
    --genomes /dir/to/FASTAs                       \
    --references /dir/to/reference/annotation(s)   \
    --emapper /dir/to/eggnog-mapper-database/data/ \
    --output ./                                    \
    --nproc 2
```

### How to cite
[Manuscript submitted]
