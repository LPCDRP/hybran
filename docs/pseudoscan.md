# pseudoscan

## Pseudo calling Criteria

1. Reference Gene is Pseudo
    -Every annotation will also be marked pseudo unless ALL of the following criteria are satisfied:
        -If the isolate gene does NOT have a reference corresponding start and stop.
        -If the isolate gene has no internal stops and ends with a valid stop codon.
        -If the isolate gene is divisible by three, while the reference gene is NOT divisible by three.
            -> assign non-pseudo

2. All other genes
    -Annotations will be marked pseudo if ANY of the following criteria are satisfied:
        -If the isolate gene has an internal stop or ends with a non-valid stop codon.
            -> assign pseudo
        -If the isolate gene is not divisible by three - invalid reading frame.
            -> assign pseudo
        -If the isolate gene doesn't have either a reference corresponding start & stop, or a strong blastp match at the user specified thresholds.
            -> assign pseudo

Each pseudo gene in the final annotation will contain a note field describing the primary reason it was assigned the 'pseudo' tag. Also, some non-pseudo genes with interesting attributes will also contain notes describing why they are NOT marked 'pseudo', despite having a certain mutation. See below for special pseudo notes assigned and the specific conditions that trigger them:

1) "Has deletion mutation(s) compared to the reference"
    -Triggered when any gene has a reference corresponding start and stop, and the reference gene is longer than the isolate gene.

2) "Has insertion mutation(s) compared to the reference"
    -Triggered when any gene has a reference corresponding start and stop, and the reference gene is shorter than the isolate gene.

3) "Has a frameshift mutation leading to a delayed stop codon"
    -Triggered when a pseudo gene has a reference corresponding start (no stop), has a valid reading frame (divisible by three) and does not contain internal/invalid stop codons.
    -This note is used to find and track gene fusion events, which is denoted by the 'geneA::geneB' nomenclature in the gene name field
