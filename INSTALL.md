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
    * Required database: bactNOG v4.5.0
    * To download: `download_eggnog_data.py -y bactNOG`

All dependencies are available within the `environment.yml` within an anaconda environment

## How to install

```
# Set up the anaconda environment with the dependencies
conda env create -n hybran -f environment.yml
conda activate hybran

# install hybran
python setup.py install
```

