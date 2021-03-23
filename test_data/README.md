# Benchmarkdata

## training data:
- Paris downloaded data from [ChiRA analysis history](https://rna.usegalaxy.eu/u/videmp/h/paris-analysis)
  - 3 replicat of mouse ES, human HEK293T celline
- Possible other are: Splash and Liga-Seq 




## What we want:
- ID 
- Interaction region
- secondary sturcture (would be good)

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



##### RRI Benchmark set:
- 158 RRIs with the position of the interaction within the mRNA
- We can add the instances that Martin found:
MicA_NC_000913 b0411_NC_000913_-200+100 C11G&G-46C [18]
RybB_NC_000913 b0805_NC_000913_-200+100 C2G&G-71C [18]
RybB_NC_000913 b2594_NC_000913_-200+100 C2G&G4C [18]
ArcZ_NC_000913 b2741_NC_000913_-200+100 C70G&G-103C [21]
ArcZ_NC_000913 b3546_NC_000913_-200+100 C69G&G-10C [22]
CyaR_NC_000913 b1824_NC_000913_-200+100 G32C&C-11G [3]
FnrS_NC_000913 b2531_NC_000913_-200+100 C47G&G-3C [3]
RyhB_NC_000913 b3365_NC_000913_-200+100 C45G&G5C [3]
MicF_NC_003197 STM0366_NC_003197_-200+100 C6G&G-31C [24]
RybB_NC_003197 STM0999_NC_003197_-200+100 C2G&G-39C [25]
RybB_NC_003197 STM1572_NC_003197_-200+100 C2G&G25C [25]
RybB_NC_003197 STM1995_NC_003197_-200+100 C2G&G19C [25]

Spot42_NC_000913 b3962_NC_000913_-200+100 G49C&C21G,U50A&A20U,A51U&U19A[30]
Spot42_NC_000913 b1302_NC_000913_-200+100 G55C&C-33G,G56A&C-34U,A57C&U-35G[31]
Spot42_NC_000913 b2715_NC_000913_-200+100 G5A&C-25U,G6C&C-26G,U7G&A-27C[31]
Spot42_NC_000913 b3224_NC_000913_-200+100 G5A&C-54U,G6C&C-55G,U7G&A-56C[31]
Spot42_NC_000913 b1901_NC_000913_-200+100 G5A&C-34U,G6C&C-35G,U7G&G-36U[32]
GcvB_NC_000913 b1130_NC_000913_-200+100 C158G&G-13C,U157A&A-14U,G156C&C-15G ,U155A&A-16U,C154G&G-17C[34]
ArcZ_NC_000913 b1892_NC_000913_-200+100 U78A&G-60U,U77A&A-59U,G76G&C-58C ,U75C&A-57G,G74C&U-56G,G73U&C-55A [35] OxyS_NC_000913 b1892_NC_000913_-200+100 A69U&U-18A,A68A&U-17U,U67C&A-16G ,A66C&U-15G,A65U&U-14A[35]
RybB_NC_000913 b0721_NC_000913_-200+100 U12A&A-21U,C13G&G-22C[36]
RyhB_NC_000913 b0721_NC_000913_-200+100 U51A&A-15U,A50U&U-14A[36]
Spot42_NC_000913 b0721_NC_000913_-200+100 G13C&C-53G,G14C&C-54G[36]
FnrS_NC_000913 b0755_NC_000913_-200+100 C47A&G-4U,U48A&A-5U,U49G&A-6C[5]
FnrS_NC_000913 b1479_NC_000913_-200+100 U57A&A-13U,U58G&A-14C,U59A&A-15U[5]
FnrS_NC_000913 b2153_NC_000913_-200+100 G5U&C-18A,G4C&C-19G[5]
FnrS_NC_000913 b2303_NC_000913_-200+100 G5U&C-6A,G4C&C-5G[5]
RyhB_NC_000913 b0592_NC_000913_-200+100 G53U&C-7A,C54U&G-8A,U55C&A-9G[38]
RyhB_NC_000913 b2155_NC_000913_-200+100 C47A&G-47U,A48U&U-48A,C49A&G-49U ,A50C&U-50G[38]
Spot42_NC_000913 b1761_NC_000913_-200+100 C46G&G86C,C48G&G88C[3]
GcvB_NC_003197 STM3930_NC_003197_-200+100 U84G&A-15C,G85A&C-16U,U86A&A-17U ,U87C&A-18G,U88A&A-19U[41]








#### Criteons to select data:
- Do we need Direct Duplex Detection (DDD) methods data to increas our trainingssets? Methods are LIGR-seq PARIS and SPLASH. Or new method  RIC-seq


