# :cherries: Cherri (Computational Help Evaluating RNA-RNA interactions) :cherries:

## Progamm Idea:
We are developing a tool that distinguishes between biological relevant or real RRIs and RRIs which most likely would not occurs in nature. Giving the tool a RRI prediction it will evaluate weather a prediction is likly to be a real prediction. 
Sofare we developed models base on the human and mouse RNA-RNA interactiom. 
Futher more the option of training a novel model basted on RNA-RNA interactom data is given.


## Instalation

Cherri developed in Linux and tested on Ubuntu (18.04 LTS). For the installation on your conda is requierd. 


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the Python 3.8 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be evocable on the command line via (possibly in a more advanced version):

```
$ conda --version
conda 4.10.1
```

### Install Cherri Conda package 


Cherri is available as a Conda package [here](https://anaconda.org/bioconda/).

We recommend to create a new Conda environment inside which we will then install cherri:

```
conda create -n cherri python=3.8 cherri -c conda-forge -c bioconda
```
Or install it into your faivorint existing enviorment
```
conda install -c bioconda cherri
```

If you are experiencing problems while running `conda install -c bioconda cherri` (e.g. complaints about conflicting dependencies), the following commands should do the trick:

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

This tells conda to explicitly look for packages in the specified channels, stored in the `.condarc` [conda configuration file](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html).


Now Cherri should be available inside the environment:


```
cherri -h
```

### Manual installation

To manually install Cherri, we first create a Conda environment. Once inside the environment, we need to install the following dependencies:

```
conda create -n cherri python=3.8 -c conda-forge -c bioconda
conda activate cherri
```

And than install the following dependencies:

```

conda install -c conda-forge markdown
conda install -c bioconda ucsc-twobitinfo
conda install -c bioconda ucsc-twobittofa
conda install -c bioconda IntaRNA

conda create -n test_biofilm -c smautner -c conda-forge -c bioconda biofilm scipy networkx

pip install eden-kernel
```

Cherri was tested with the following versions: . 
Finally, to install the tool itself, we simply clone the repository and execute the installation script inside the folder:

```
git clone https://github.com/teresa-m/Cherri.git
cd Cherri
python -m pip install . --ignore-installed --no-deps -vv 
```

Now we can run Cherri from any given folder (just remember to re-activate the environment once you open a new shell):

```
cherri -h
```



### Usage
You can use Cherri in two modes. The first mode to prdict weather your RRI prediction biolocical relevant. You can used the current model for human or mouse data. If you would like build a model based on novel RRI interactome data you can use the seciond mode.

#### First mode: evaluation of RRIs
Based on a tabluar file containing chromosomal position data of the RRIs, it classify if the interaction region is likly to be a relevant one.

#### Second mode: build new Cherri model
Using Chira RRI output data the build model mode will generate a prediction classifier. 

### Backround of scripts which are wrapped
The following scripts are needed to build the trainings data to build a model.


#### First mode: evaluate_instance.py
Input would need the postion of the interaction. And the genome file and a file stroing the chromosome length is needed to rund a evaluation. The evaluation script calls one after a other differen 


#### Input format
The instances which should be tested should be in a tabluar format. You should specify a header lines, with the chromosome interaction start interaction end and strand of the two interacting partners. 
```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```
With a specific setting only 'positve' instances are computed. If no occupyed regions are not given by using none this can be specifyed. If the occpyed regions for human or the mouse model used for training should be used specify mouse or human. 
#### example call
```
python evaluate_instance.py -i1 /vol/scratch/data/RRIs/test/test_evalueat_rris.cvs -i2 none -g /vol/scratch/data/genomes/hg38_UCSC_20210318.2bit  -o /vol/scratch/data/test/ -c 150 -n test_data -l /vol/scratch/data/genomes/hg38_Info.tab
```

#### Input Parameter
- i1 | RRIs_table: table containg all rris that should be evalutated in the corrct fromat
- i2 | occupyed_regions: path to occupyed regions file or specify: human, mouse or none 
- g | genome_file: path to 2bit genome file
- o | out_path: path to folder all output folder of each step of the data preparation
- c | context: how much context should be added at left an right of the sequence
- n | experiment_name: name of the datasoruce of RRIs
- l | chrom_len_file: tabular file containing chrom name \t chrom lenght for each chromosome
- p | param_file: IntaRNA parameter file

#### Output 
Within the output path a folder is created with the date of today used as name. 
Within this folder the IntaRNA call results are stored in the positive_instance folder and the features are stored in the feature_files folder.


#### Second mode: create_trainings_data.py
Is a wrapper script that will take RNA-RNA interacting subsequences data generated by CheRRI. The wrapper will call the following scripts afte each other: 
1. find_occupied_regions.py: outputs in the folder 'occ_out' a trusted rri file and a occpyed_region object.
2. generate_pos_neg_with_context.py: outputs in the folder 'pos_neg_data/' the positive and negative instances 
3. get_features.py: outputs in the 'feature_files/' folder the postive and negative tablular feature files


#### example call
```
python create_trainings_data.py -i1 /vol/scratch/RNA_RNA_binding_evaluation/test_data/training/Paris/ -i2 /vol/scratch/data/ -g /vol/scratch/data/genomes/hg38_UCSC_20210318.2bit -r test_rep1.tabular test_rep2.tabular test_rep3.tabular  -o /vol/scratch/data/test/ -c 150 -n test_data -l /vol/scratch/data/genomes/hg38_Info.tab
```

#### Input Parameter
- i1 | RRI_path: path to folder storing the ChiRA interaction summary files
- g | genome_file: path to 2bit genome file
- r | list_of_replicats: list ChiRA interaction summary files names of all replicats
- o | out_path: path to folder all output folder of each step of the data preparation
- c | context: how much context should be added at left an right of the sequence
- n | experiment_name: name of the datasoruce of RRIs
- l | chrom_len_file: tabular file containing chrom name \t chrom lenght for each chromosome
- p | param_file: IntaRNA parameter file

#### Output 
- All data sources are in the 'Cherri_model_build/' folder. 





### IntaRNA params for call

To generate the current features we need:

| param  | value  | decription | 
|---|---|---|
| outMode |  C  | Output style (C=tabluar format) |
| seedBP |  5  | the number of base pairs within the seed |
| seedMinPu |  0  | the minimal unpaired probability of each seed region in query and target |
| accW |  150  | sliding window length (0=global folding) |
| acc |  N/C  |  To globally turn off accessibility consideration: turn off/on |
| outMaxE |  -5  | maximal energy for any interaction reported |
| outOverlap |  B  | overlapping of interaction sites of Suboptimal allowed (B:both) |
| outNumber |  5  | generate up to N interactions for each query-target pair |
| seedT/QRange |  positive interaction |  |
| q/tAccConstr |  negative interaction  |  |
| intLenMax |  50 | restrict the overall length an interaction |
| temperature |  37  |  |
| intLoopMax |  3  | number of unpaired bases between intermolecular base pairs |



### Interaction input

To generate the current features we need:

| description  | coloum name  | 
|---|---|
| start positon of RRI target  | 'start1'  |
| end positon of RRI target   |  'end1' | 
|  start positon of RRI query  |  'start2' |
| end positon of RRI query   |  'end2' | 
| dotbrackad of interaction  | hybridDP  |
| RRI sequence target&query  |  subseqDP | 
| MFE  |  E | 
|  seeds starts listed for target |  seedStart1 | 
| interacting target sequence  |  target | 
|  interacting query sequence | query |


## Script explanations:

### find_trusted_RRI.py
Here we search for trusted RRIs, so RRIs which can be found in all replicates. In a first filter step only uniquely mapped RRIs are taken. Than RRI sequence partners in all replicas are found, using a overlap threshold. Output are the Chira input tables, now containing only the trusted RRIs. Only one of the sequence RRI pairs is added to the output. 

#### example call
```
python find_trusted_RRI.py -i /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/data/data_Chira/training/Paris -r test_rep1.tabular test_rep2.tabular test_rep3.tabular -o 0.6 -d /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/ -n test_paris
```

#### Input Parameter
- input_path: path to folder storing all input data
- list_of_replicats: list_of_replicates
- overlap_th: path output reposetory
- experiment_name: name of the data soruce of positve trusted RRIs

#### Output 
- trusted RRIs in tabulat format
- pickled list of avg overlap percent
- pickled list of avg overlap lenght


### plot_avg_overlaps.py 
The scirpt takes the lists of avg overlap percent and overlap lenght form the find_trusted_RRI.py script. Thie Idea ist that for differen overlap thresholds will be plotted into a box plot. At the moment for overlap thresholds ['0.3', '0.4', '0.5', '0.6', '0.7']. Therefor the find_trusted_RRI.py should be called with this overlap thresholds.

#### example call
```
python plot_avg_overlaps.py -i /vol/scratch/data/plots/overlap_distibution/ -n rri_overlap_plots
```

#### Input Parameter
- input_path: path to folder where input files
- experiment_name: name of the data soruce of positve trusted RRIs


#### Output 
- two boxplots



### find_occupied_regions.py
Given the RRI information tables form Chira and RNA-Protein binding positions a Interlab object is build. The occuped information can be used to mask parts of the genome and therefore enale to select negative interaction regions. This reagions are not part of interaction in nature. 

#### example call
```
python find_occupied_regions.py -i1 /vol/scratch/data/RRIs/Paris/ -r test_rep1.tabular test_rep2.tabular test_rep3.tabular -o /vol/scratch/data/RRIs/
```

#### Input Parameter
- RRI_path: path to folder storing all RRI data (tabular)
- rbp_path: path to RBP side data file (bed format)
- list_of_replicats: list having filenames of all replicats
- out_path: path to folder storing outputfiles
- overlap_th: overlap threshold
- score_th: score threshold

#### Output 
- outputs a pickels Interlab object and prints the path to the file
- since it calles trusterd RRIs script the output will also be a file for trusted RRIs once with a lower score threshold to bild the occupyed regions from and a score filterd file, which can be used to build the build postive and negative prediction instances. 



### generate_pos_neg_with_context.pl
Given trusted RRI and the occupyed regions a given context is appended. Than positive sequences are computed by calling IntaRNA specifiying the the seed region with the trusted RRI regions. The negative interactions are computed by calling IntaRNA given regions, which should not be in the interaction 'occupied regions'. 

**sofare we can only coumpute the this data for HEK293T RRIs**

#### example call
```
python generate_pos_neg_with_context.py -i1 /vol/scratch/data/trusted_RRIs/test_paris_overlap_0.6.cvs -i2 /vol/scratch/data/RRIs/20210809_occ_out/occupied_regions.obj -d /vol/scratch/data/pos_neg_data_context/ -g /vol/scratch/data/genomes/hg38_UCSC_20210318.2bit -n test_paris_HEK293T  -c 10
```

#### Input Parameter
- input_rris: path to file storing all positve trusted RRIs
- input_occupyed: path to file storing occupied sides (InterLab obj)
- output_path: path output reposetory
- experiment_name: name of the datasoruce of positve trusted RRIs
- genome_file: path to 2bit genome file
- context: how much context should be added at left an right of the sequence
- pos_occ/ no_pos_occ: if postived instances should have the occupyed positions blocked with the IntaRNA call
- block_ends: number of blocked nucleotides at the start and end of the sequence
- no_sub_opt: # of ineractions IntraRNA will give is possible
- chrom_len_file: tabular file containing chrom name \t chrom lenght for each chromosome
- param_file: IntaRNA parameter file

#### Output 
- Positive and negative datasets stored in tabular format.



### get_negative_dataset.py
Generate the positive and negative dataset for the trusted RRIs. A context around the interaction side can be specified and what kind of shuffling should be used to generate the negative data. The negative instance is chosen form a random number of options, which has the closed energy to the positive RRIs. 

#### example call
```
 python get_negative_dataset.py -i /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/test_paris_HEK293T_overlap_0.6.cvs -d /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/ -g /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/data/genomes/hg38_UCSC_20210318.2bit -n test_paris_HEK293T -k 3 -s 5 -cm together -c 10
```

#### Input Parameter
- input_file: path to file storing all positve trusted RRIs
- output_path: path output reposetory
- experiment_name: name of the datasoruce of positve trusted RRIs
- kind_of_shuffel: seqence mononucleotide (1) or sequence denucleotide (2) or bp mononucleotide (3) shuffling
- shuffle_no_seq: how often is the positive sequence shuffled
- context_method: select the context method  if context should not be added (non), if it should be shuffled sepatatly (separat), or together (together) with the sequence
- context: how much context should be added at left an right of the sequence

#### Output 
- Positive and negative datasets stored in tabular format. The table contains the all informatio of the trused RRI instance and the results of the IntaRNA call. 


### plot_tRRIs.py
For the postive and negative datasets overview plots are generated. 

#### example call
```
python plot_tRRIs.py -i1 test_paris_HEK293T_context_method_together_shuffling_method_3_pos_RRI_dataset.csv -i2 test_paris_HEK293T_context_method_together_shuffling_method_3_neg_RRI_dataset.csv -i3 /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/ -o /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/
```

#### Input Parameter
- input_pos: file name of positive dataset in tabular format (, separted)
- input_neg: file name of negative dataset in tabular format (, separted)
- save_path: directory where the plots will be stored in a plot folder

#### Output 
- RNA historam plot
- Energy distibution plot


### get_features.py
Here for a given input and a given reature set the features, this features are stored in a tabular format. 


#### Our List of features
- E : Minimum free energy (overall interaction energy)
- E_hybrid : energy of hybridization only = E - ED1 - ED2
- maxED: maximal energy of ED1 and ED2

- no_bps : number of base pairs within the interaction
- GC_content: GC content of the interaction side
- max_inter_len: the maximum interaction side calculated by the length of the target sequence and query sequence length
- no_seeds: numbers of possible seeds within the interaction

- inter_len_normby_bp: the maximum interaction length divided by the number of base pairs with the interaction 
- bp_normby_inter_len: number of base pairs within the interaction divided by the maximum length at the interaction
- mfe_normby_GC: MFE divided by the GC content
- max_ED_normby_GC: max_ED divided by the GC content
- E_hybrid_normby_GC: E_hybrid divided by the GC content
- complex_target: sheenon entropy of target sequence
- complex_query: sheenon entropy of query sequence

#### example call
```
python get_features.py -i ../output/paris_HEK293T_06_context_method_together_shuffling_method_2_pos_RRI_dataset.csv -f E no_bps GC_content mfe_normby_GC -o /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/input_features/2
```

#### Input Parameter
- input: path to input file
- feature_set_list: list of all featurs that should be summarized in the output
- output_file: file path where the output table should be stored

#### Output 
- tabular file having all features given via the feature_set_list



### wrapper_feature_generation.py
Building a table and a call.sh script for a combination of different negative data sets and feature combinations. This script helps to automatize a grid search to to find the best data and parameters for a ML model.

#### example call
```
python wrapper_feature_generation.py -i /vol/scratch/data/pos_neg_data -s human_mouse -f 'E no_bps GC_content max_inter_len inter_len_normby_bp bp_normby_inter_len mfe_normby_GC no_seeds complex_target complex_query' 'no_bps GC_content max_inter_len' 'bp_normby_inter_len mfe_normby_GC complex_target complex_query' -o /vol/scratch/data/feature_input_human_mouse
```

#### feature set combiations:
- all: 'E no_bps GC_content max_inter_len inter_len_normby_bp bp_normby_inter_len mfe_normby_GC no_seeds complex_target complex_query'
- small easy: 'no_bps GC_content max_inter_len'
- small complex 'bp_normby_inter_len mfe_normby_GC complex_target complex_query'

#### Input Parameter
- input_dir: path to input files, the file names are currently stored within the script
- feature_set_list: list of all featurs combinations, each list should be given in ''
- input_set: prefix indicating the experiment and subest to be analyzed
- output_dir: dir where feature tables will be stored

#### Output 
- tabular file and call.sh file. The tabel will give a overview for the possible combitoantion of feature sets. Once the call.sh file is perfromed, all of the modle input tabels or vectors are generated. They can be used to train different models.




### run_different_feature_combinations.py
calls the training.py script for different for all negative dataset and feature combinations. It takes the compintions for the overview tabular. 

#### example call
```
python run_different_feature_combinations.py -i /vol/scratch/data/feature_input_snRNA/ -o /vol/scratch/output/modle_test/ -e snRNA
```

#### Input Parameter
- input_dir: "path to input files
- out_dir: path to output dir
- experiment: prefix indicating the experiment and subest to be analyzed

#### Output 
- tabular file contiaing the feature and experiment inforamiton and the AUC and STD of the different models.



### plot_heatmap.py
plots a heat map for AUCs of feature set and negative data combinations. 

#### example call
```
python plot_heatmap.py -i /vol/scratch/output/modle_test/lncRNA_modle_AUC.csv -o /vol/scratch/output/modle_test/ -c 10 -d paris_lncRNA
```

#### Input Parameter
- input_file: path to input files
- out_dir: path to output dir
- context: context
- data_name: name of dataset

#### Output 
- heat map and updated table.



### python training.py


#### example call
```
python training.py -i /vol/scratch/data/feature_input_human_mouse//1a  -n /vol/scratch/data/feature_input_human_mouse//1b  -d /vol/scratch/output/modle_test/
```

#### Input Parameter
- i: positive feature set
- n: negative feature set
- f: output dir where feature tables will be stored

#### Output 





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


#### Sebastians feature selection
- MFE
- Maximal length of the two interacting subsequence
- Number of base pairs within the top 1 RRI 
- Maximal length of an interacting subsequence normalized by the number of base pairs within the top 1 RRI
- Number of base pairs within the interaction vs. the normalizedmaximal length of the top 1 RRI
- Maximal energy ED to make one of the interacting subsequencesof the top 1 RRI 
- GC-content within interaction side
- Requencies of the minimum free energy normalizedby the GC-content
- Number of possible seeds


