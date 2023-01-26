# Comparing annotation results from different versions of hybran

Two helper scripts in the `tests` folder make this possible.
Use `tests/gbk2coords` to create a tabular version of the final genbank output, then use `tests/compare-summaries` to compare two of them.
The latter script is a small wrapper around `diff` that removes the numbers from ORF#### gene names, as these numbers can vary and appear as spurious differences.

Example:

```
./tests/gbk2coords ../verA-results/sample.gbk > ../verA-results/coords.txt
./tests/gbk2coords ../verB-results/sample.gbk > ../verB-results/coords.txt
./tests/compare-summaries ../ver{A,B}-results/coords.txt | less
```
