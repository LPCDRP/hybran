# AnnoTUB

Annotation pipeline that provides genomic annotation 
for *Mycobacterium tuberculosis* complete genomes

### Dependencies
* [RATT](http://ratt.sourceforge.net/)
* [Prokka](https://github.com/tseemann/prokka)
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [python2.7](https://www.python.org/downloads/release/python-2715/)
* [Biopython](https://biopython.org/wiki/Download)
* [numpy](https://sourceforge.net/projects/numpy/files/NumPy/) (is also available on `pip` and [anaconda](https://anaconda.org/anaconda/numpy))
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [MCL](https://github.com/JohannesBuchner/mcl)
* [NCBI-BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Required files
* *M. tuberculosis* *de novo* assembled, complete genome
* [H37Rv CDS FASTA](https://www.ncbi.nlm.nih.gov/nuccore/AL123456.3)

### How to run
**RATT_HOME variable must be defined and exported!**
`make -f ratt_prokka.mk GENOMES=dir/to/fastas H37Rv_FASTA=H37Rv-CDS.fasta`

`make -f clustering.mk annotation_dir=dir/to/gffs`

