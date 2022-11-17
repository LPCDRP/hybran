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

All dependencies are available within the `environment.yml` within a [conda](https://docs.conda.io/en/latest/miniconda.html#installing) environment

## How to install

The easiest way to install `hybran` is through the conda package manager:

```
conda install -c bioconda hybran
```

### Installation from Source / Development Installation

* Download and extract the [latest release](https://gitlab.com/LPCDRP/hybran/-/releases/permalink/latest).

* Install the dependencies using [miniconda](https://docs.conda.io/en/latest/miniconda.html#installing) (or manually)
```
# From the hybran directory you downloaded, set up the conda environment with the dependencies
conda env create -n hybran-dev -f environment.yml
conda activate hybran-dev
```

* Install hybran
```
pip install --editable .
```

If you run into any trouble, please let us know over at our [issues page](https://gitlab.com/LPCDRP/hybran/-/issues).
