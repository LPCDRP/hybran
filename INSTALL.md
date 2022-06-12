# Installation

## Dependencies
* [RATT](http://ratt.sourceforge.net/)
* [Prokka](https://github.com/tseemann/prokka)
* [EMBOSS (seqret)](http://emboss.sourceforge.net/download/)
* [Biopython](https://biopython.org/wiki/Download)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [MCL](https://github.com/JohannesBuchner/mcl)
* [NCBI-BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)

All dependencies are available within the `environment.yml` within an anaconda environment

## How to install

* Download and extract the [latest release](https://gitlab.com/LPCDRP/hybran/-/releases/permalink/latest).

* Install the dependencies using [miniconda](https://docs.conda.io/en/latest/miniconda.html#installing) (or manually)
```
# From the hybran directory you downloaded, set up the conda environment with the dependencies
conda env create -n hybran -f environment.yml
conda activate hybran
```

* Install hybran
```
python setup.py install
```

