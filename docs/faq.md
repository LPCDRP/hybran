
# Special Considerations

## Reannotating your reference genome with Hybran

You will find out that *ab initio* predictions are made in your sample that would have also been found in your reference genome.
To get a clearer view of the genes present only in your sample, you may want to consider reannotating your reference with Hybran.
In this case, you'll find it helpful to use a different ORF prefix with which to name unknown genes, as well as RATT transfer type "Assembly".

:warning: If you want to maintain the original locus tag prefix, Your reference genome fasta file **must** be named using it.
Additional genes will be assigned locus tags counting from the highest number observed in the annotation.
The advantage is that the original locus tags will be maintained, but the disadvantage is that you'll introduce unapproved new labels for the additional genes you find.

```
ln -s H37Rv.fasta Rv.fasta # This strain's locus tag prefix is Rv
hybran --genomes Rv.fasta --references H37Rv.gbk --dedupe-references --orf-prefix RvMTB --ratt-transfer-type Assembly ...
```

## Annotating multiple genomes

Hybran assigns numbered generic names to unknown genes (HYBRAN1234, for example) so that you can match gene names across samples.
If you have multiple genomes to annotate and want to take advantage of this, the best approach is to annotate them simultaneously-- Hybran can accept multiple fasta files as input.
In the example below, all the genomes in the `assemblies` folder will be annotated.

```
hybran --genomes assemblies -r H37Rv.gbk -o annotations ...
```

If you have annotated several samples already and later receive a new sample to annotate, you can still do so without reannotating everything else, and all previously observed `HYBRAN####` names will be used consistently.
Just pass all your previously annotated genomes as reference genomes together with your primary reference.
But make sure to set your primary reference as the `--first-reference` since it's the only one used for initial reference matching to the *ab initio* calls.

So if your annotations are in the folder `annotations`, as created above, you can create a file of file names (fofn) for your references and use it as follows:

```
ls $PWD/annotations/*.gbk $PWD/H37Rv.gbk > refs.fofn
hybran --genomes assemblies/new-sample.fasta --references refs.fofn ...
```

## Large Number of Reference Annotations

If you are using a large number of references, you may get an error during the annotation transfer step like this:

```
ERROR: mummer and/or mgaps returned non-zero
```

If so, you should recompile mummer using `make CPPFLAGS="-O3 -DSIXTYFOURBITS"` as described at <https://sourceforge.net/p/mummer/mailman/message/34069398/>.
