# Benchmarking CheRRI

We compiled a [data set of of experimentally verified RNA-RNA interactions for humans](benchmark-data.xlsx).
In the following, we list the publications and the respective ids within the benchmark data set.

| ids | doi of publication |
| --- | ---- |
| um.1-38 | [Umu & Gardener, 2017](https://doi.org/10.1093/bioinformatics/btw728) |

Each RNA pair comes with support from experiments where the interaction is formed.

To benchmark CheRRI, we did the following for each pair of RNAs:

1. run prediction tool XYZ to identify a putative interaction site (see [benchmark-calls.txt](benchmark-calls.txt))
2. extract the interaction site: if the site shows an overlap of at least 5 nt with the support information rate it a "correct prediction", otherwise "wrong". this defines the *gold standard*. (see [benchmark-stats.R](benchmark-stats.R))
3. run CheRRI on the extracted interaction site: compare CheRRI's rating with the gold standard to assess whether CheRRI's prediction is correct (TRUE POS|NEG) or wrong (FALSE POS|NEG).

We did this for the following tools

- IntaRNA version 3.3.2 (2022-09-13) [intaRNA-3.3.2.tar.gz](intaRNA-3.3.2.tar.gz)
- RIblast version 1.2.0 (2019-11-02) [RIblast-1.2.0.zip](RIblast-1.2.0.zip)
- RIsearch2 version 1.2 (2021-06-16) [RIsearch-1.2.tar.gz](RIsearch-1.2.tar.gz)

The interaction predictions, their correctness assignments and resp. CheRRI evaluations are provided in [overlap.csv](overlap.csv).

Based on this, we mapped the interacting regions to respective genomic positions (hg38) to generate the [cherri_input.csv](cherri_input.csv) file, which holds both each region of interest as well as its correctness and source annotation.

Given this, we provide some statistics on the correctness of CheRRI's predictions.
