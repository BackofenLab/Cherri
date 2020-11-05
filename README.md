# RNA_RNA_binding_evaluation

## Project Idea:
Predict if a given RNA-RNA interaction could be a real biolocal one. Generate a ML-Model that can do this prediction. For this we will have the following work plane

### Dataset retreaval
Options for the dataset retreaval
1. intraRNA benchmark: It has 149 validated sRNA-mRNA interactions. Find in in the [IntaRNA-benchmark](https://github.com/BackofenLab/IntaRNA-benchmark). The varified interaction.tsv holds lists the sRNA-mRNA including there sorce publication and if known compensatory mutaions. 
2. RNA Interaction Database [RNAInter](http://www.rna-society.org/rnainter/). It contains RNA-RNA interaction data only for homo sapiens for stong interaction and 2254 week experimantaly verfiyed interactions. RRI included for the follwing dbs: RISE, LncRNA2Target v2.0, VIRmiRNA, LncACTdb 2.0, NPInter v3.0, OncomiRDB, ncRDeathDB, miR2Disease, sRNATarBase, MNDR v2.0, LncRNADisease 2.0, VmiReg. Thea came to have 6 007974 RRIs. (Table 1. Overview of curated interaction data from 35 resources)
3. [RISE](http://rise.life.tsinghua.edu.cn/downloads.html) databas contains RRI from human, mouse, yeast, E.coli, S.enteria, or 10 different cell lines. The data is stored in BEDPE file format. (-> data included in RNAInter) 
4. [miRTarBase](http://mirtarbase.cuhk.edu.cn/php/index.php) has data itrcgated form miRBase, NCBI Entrez gene, NCBI RefSeq, SomamiR, HMDD, TransMir, miRSponge, CMEP, Gene Expression Omnibus (GEO) (32), The Cancer Genome Atlas (TCGA). miRTarBase contains 479,340 curated MTIs between 4,312 miRNAs and 23,426 target genes. (-> data included in RNAInter, but maybe a older version!) 
5. [sRNATarBase 3.0](http://ccb1.bmi.ac.cn:81/srnatarbase/) has 475 validated sRNA-mRNA interactions. Is not working stabaly... (-> data included in RNAInter) 
6. Studys:
- RNA-RNA interactions in E. coli using Clash https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3725-3. maybe intersting for later
- https://academic.oup.com/bioinformatics/article/33/7/988/2557685
- https://onlinelibrary.wiley.com/doi/full/10.1002/wrna.1600
- https://academic.oup.com/bioinformatics/article/35/16/2862/5262219
- https://academic.oup.com/bioinformatics/article/35/17/3199/5298307
- https://www.sciencedirect.com/science/article/pii/S109727651630106X
- https://www.sciencedirect.com/science/article/pii/S0092867418309486


#### RNAInter
1. **Downlad** All RNA-RNA interactions gives a file with follwoing colums:
- RNAInter ID: unique identifier for each entry in RNAInter database
- Interactor1: interactor1 in current entry
- ID1: ID of interactor1
- Category1: category of interactor1
- Species1: organism name of interactor1
- Interactor2: interactor2 in current entry
- ID2: ID of interactor2
- Category2: category of interactor2
- Species2: organism name of interactor2
- Score: the integrative confidence score of the current entry
2. **Browse** folders RRI sorted by RNA e.g. mRNA, miRNA.. or species. Files have the follwoing colums:
- RNAInter_ID	
- Interactor1	
- Category1	
- Species1	
- Interactor2	
- Category2	
- Species2	
- Score
3. **Search** can do a search based on a keyword Categry interactiontype ... Problem is that we need a keyword!
The same file columes, but in the online version you have the coulms Details. 

RIscoper- tool for RRI extraction from literature (RNAInter)




#### Criteons to select data:
- Do we need Direct Duplex Detection (DDD) methods data to increas our trainingssets? Methods are LIGR-seq PARIS and SPLASH. Or new method  RIC-seq




### select features

- nucleotied-contens, nucleodited-skews and sequence complexety measurs
- featurs for lncRNAs: more cis than trans binding (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0150353)


Feature from RNAz:
- G+C content, 
- the A/(A+U) ratio, 
- the C/(C+G) ratio, 
- all 16 dinucleotide frequencies
- the length of the sequence scaled to the interval [0,1]

- regression for estimating the mean free energy was trained to learn energy per nucleotide
- standard grid search approach was used to find optimal combinations for SVM parameters.



### select ML-Method
