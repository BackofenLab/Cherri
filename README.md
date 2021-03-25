# RNA_RNA_binding_evaluation

## Project Idea:
Predict if a given RNA-RNA interaction could be a real biolocal one. Generate a ML-Model that can do this prediction. For this we will have the following work plane






## select features

- nucleotied-contens, nucleodited-skews and sequence complexety measurs
- featurs for lncRNAs: more cis than trans binding (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0150353)


#### Features according to the Antrag: 
- composition of the predicted duplex (such as interaction energy, number of
intramolecular base pairs, number of unpaired bases, maximum loop length, seed length,
average/maximum position-wise probability of an interaction etc). 
We will also use
- information about absolute and relative (compared to the context) conservation of the
associated binding sites as well as the p-value to get a hit in the genome that has a higher
energy (i.e., against the background of possible interactions for the considered RNA) 

#### Feature from RNAz:
- G+C content, 
- the A/(A+U) ratio, 
- the C/(C+G) ratio, 
- all 16 dinucleotide frequencies
- the length of the sequence scaled to the interval [0,1]

- regression for estimating the mean free energy was trained to learn energy per nucleotide
- standard grid search approach was used to find optimal combinations for SVM parameters.


#### Sebastians Feature selection
- MFE
- maximal length of the two interacting subsequence
- number of base pairs within the top 1 RRI as predicte
- maximal length of an interacting subsequence normalized by thenumber of base pairs within the top 1 RRI
- number of base pairs within the interaction vs. the normalizedmaximal length of the top 1 RRI
- maximal energyEDto make one of the interacting subsequencesof the top 1 RRI 
-  GC-content within the top 1 RRI as predicte
- minimum free energy and GC-content of thetop 1 RRI 
- requencies of the minimum free energy normalizedby the GC-content
- number of possible seeds within the top 1 RRI




### select ML-Method
