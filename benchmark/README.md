# Benchmarking CheRRI

We compiled a data set of of experimentally verified RNA-RNA interactions for humans.
In the following, we list the publications and the respective ids within the benchmark data set.

| ids | doi of publication |
| --- | ---- |
| um.1-38 | [Umu & Gardener, 2017](https://doi.org/10.1093/bioinformatics/btw728) |

Each RNA pair comes with support from experiments where the interaction is formed.

To benchmark CheRRI, we did the following for each pair of RNAs:

1. run prediction tool XYZ to identify a putative interaction site
2. extract the interaction site: if the site shows an overlap of at least 5 nt with the support information rate it a "correct prediction", otherwise "wrong". this defines the *gold standard*.
3. run CheRRI on the extracted interaction site: compare CheRRI's rating with the gold standard to assess whether CheRRI's prediction is correct (TRUE POS|NEG) or wrong (FALSE POS|NEG).

Given this, we provide some statistics on the correctness of CheRRI's predictions.
