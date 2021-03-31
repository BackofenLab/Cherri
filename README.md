# RNA_RNA_binding_evaluation

## Project Idea:
Predict if a given RNA-RNA interaction could be a real biolocal one. Generate a ML-Model that can do this prediction. For this we will have the following work plane



## select features

#### Feature selecton tools
- nucleotied-contens, nucleodited-skews and sequence complexety measurs
- featurs for lncRNAs: more cis than trans binding (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0150353)
- feature selection using [Relief](https://doi.org/10.1016/j.jbi.2018.07.014)


#### Features according to the Proposal: 
- composition of the predicted duplex:
1. interaction energy
2. number of intramolecular base pairs, 
3. number of unpaired bases, 
4. maximum loop length, 
5. seed length,
6. average/maximum position-wise probability of an interaction etc
We will also use
- information about absolute and relative (compared to the context) conservation of the
associated binding sites as well as the p-value to get a hit in the genome that has a higher
energy (i.e., against the background of possible interactions for the considered RNA) 

#### Feature from RNAz:
compute z-scores for folding energies (shuffling and folding is too costly)
Building of the background by controling: 
- G+C content, 
- the A/(A+U) ratio, 
- the C/(C+G) ratio, 
- all 16 dinucleotide frequencies
- the length of the sequence scaled to the interval [0,1]

- regression for estimating the mean free energy was trained to learn energy per nucleotide
- standard grid search approach was used to find optimal combinations for SVM parameters.

- Features to detect structured noncoding RNAs
1. average minimum free energy z-score estimated from a dinucleotide shuffled background
2. structure conservation index or SCI (consensus folding free energy / average of the folding free energies of the single sequences)
3. normalized Shannon entropy H of the alignment as a measure for the content of evolutionary information


#### Sebastians Feature selection
- MFE
- Maximal length of the two interacting subsequence
- Number of base pairs within the top 1 RRI 
- Maximal length of an interacting subsequence normalized by the number of base pairs within the top 1 RRI
- Number of base pairs within the interaction vs. the normalizedmaximal length of the top 1 RRI
- Maximal energy ED to make one of the interacting subsequencesof the top 1 RRI 
- GC-content within interaction side
- Minimum free energy vs. GC-content 
- Requencies of the minimum free energy normalizedby the GC-content
- Number of possible seeds




### select ML-Method
