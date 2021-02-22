# Hybran

Genome annotation pipeline for *Mycobacterium tuberculosis* de novo assembled genomes

### System Requirements
* Linux OS required
* Tested on Debian 8 (jessie)
* Install time ~5 minutes

### Dependencies
* [RATT](http://ratt.sourceforge.net/)
* [Prokka](https://github.com/tseemann/prokka)
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [python2.7](https://www.python.org/downloads/release/python-2715/)
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

# install the Python 3 version of eggnog-mapper into the conda environment
pip install eggnog-mapper

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

### How to cite
[Manuscript submitted]
