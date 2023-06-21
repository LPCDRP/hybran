# Unit tests

Unit tests can be run using `pytest` in this directory.
Some data files are required to run the full suite, which will be described at some point.

# Full runs

Again, this requires some external data files.

## Options

HYBRANFLAGS
: Extra flags to pass to `hybran`

SUBMIT
: if defined, will perform the hybran runs in a batch cluster job.

RESDIR
: Resources directory where the required data files can be found. (Default: `${GROUPHOME}/data/depot/hybran/resources`)

OUTDIR
: Output directory. (Default: mm_dd (date))

NPROC
: Number of parallel threads. This applies both to local and batch runs.


## Examples

Run all unit tests on the cluster (using `qsub`) using 16 cores per run:
```
make -j SUBMIT=1 NPROC=16
```

Run a specific sample using 12 cores:

```
make 1-0006 NPROC=12
```

