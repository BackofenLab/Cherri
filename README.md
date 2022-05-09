# :cherries: Cherri (Computational Help Evaluating RNA-RNA interactions) :cherries:

## Progamm Idea:
Cherri is a tool that distinguishes between biological relevant or real RRIs and RRIs which most likely would not occurs in nature. Giving the tool a RRI prediction it will evaluate weather a prediction is likely to be a real prediction. 
Software we developed models base on the human and mouse RNA-RNA interactome. 
Further more the option of training a novel model basted on RNA-RNA interactome data is given.

[image](./plots/Cherri_workflow.pdf)


## Instalation

Cherri developed in Linux and tested on Ubuntu (18.04 LTS). For the installation on your conda is required. 


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the Python 3.8 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be accessed on the command line via (possibly in a more advanced version):

```
$ conda --version
conda 4.10.3
```

### Install Cherri Conda package **not woring jet**

Cherri is available as a Conda package [here](https://anaconda.org/bioconda/).

We recommend to create a new Conda environment inside which we will then install cherri:

```
conda create -n cherri python=3.8 cherri -c conda-forge -c bioconda
```
Or install it into your favorite existing environment
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

To manually install Cherri, you first create a Conda environment. Inside the environment, you need to install the following dependencies:

```
conda create -n cherri python=3.8 -c conda-forge -c bioconda
conda activate cherri
```

And than install the following dependencies:

```
conda install -c conda-forge scikit-learn
conda install -c conda-forge networkx
conda install -c bioconda ucsc-twobittofa
conda install -c bioconda interlap
conda install -c bioconda IntaRNA
conda install -c bioconda numpy
conda install -c bioconda pandas
conda install -c smautner eden-kernel
conda install -c smautner biofilm=0.0.88
conda install -c conda-forge python-wget
```
Or create the environment with all dependencies at ones:

```
conda create -n cherri -c conda-forge -c bioconda -c smautner scikit-learn networkx ucsc-twobittofa interlap pandas IntaRNA python-wget eden-kernel biofilm=0.0.88 python-wget

```

Cherri was tested with the following versions: XX. 
Finally, to install the tool itself, we simply clone the repository and execute the installation script inside the folder:

```
git clone https://github.com/teresa-m/Cherri.git
cd Cherri
python -m pip install . --ignore-installed --no-deps -vv 
```

You additionally need to fix the python hash seed:

```
conda env config vars set PYTHONHASHSEED=31337
```
or just for your current session!
```
export PYTHONHASHSEED=31337
```

Now we can run Cherri from any given folder (just remember to re-activate the environment once you open a new shell):

```
cherri -h
```



## Usage
You can use Cherri in two modes. The first mode to predict weather the positon of your RRI prediction is biological relevant. Precomputed models exist for human or mouse data and can be downloaded. If you would like build a model based on novel RRI interactome data you can use the second mode.

### First mode **eval** :  evaluation of RRIs
Based on a tabular file containing chromosomal position data of the RRIs, Cherri classify if the interaction region is likely to be a biolocal relevant one.

Cherris models are trained for human and mouse, which can be directly used for RRIs position evaluation. If there exists a RNA-RNA-interactome dataset for your preferred organism we recommend to train your own organism specific model using cherris train mode. This model can be than use in the eval mode for the classification of your predicted RRIs positions.

The pretrained models can be downloaded from [zenodo](https://zenodo.org/)


#### Input format of the RRIs to evaluate
The to evaluate instances should in a tabular format. You should specify a header lines, with the chromosome interaction start interaction end and strand of the two interacting partners:
```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```
With a specific setting only 'positive' instances are computed. If no additional occupied regions are specified only the once of the given input interactions are used. However, we recommend to specify the occupied regions of the trained model, to be consisted within the feature generation.


#### Input Parameter
##### required
- i1 | RRIs_table: Table containing all RRIs that should be evaluated in the correct format
- g | genome_file: Path to 2bit genome file, or used the build in download if you want the human or moues genomes
- o | out_path: Path to folder all output folder of each step of the data preparation
- l | chrom_len_file: Tabular file containing data in two columns one row for each chromosome: 'chrom name' \t 'chrom lenght'. You can speciy 'human' or 'mouse' directly.
- m | model_file: Set path to the model which should be used for evaluation
- mp | model_params: Set path to feature file of the given model 
##### optional
- i2 | occupyed_regions: Path to occupied regions object file. 
- c | context: How much context should be added at left an right of the sequence
- n | experiment_name: Name of the data source of RRIs
- p | param_file: IntaRNA parameter file. default: file in rrieval/IntaRNA\_param
- st | use_structure: Set 'off' if you want to disable structure, default 'on'
- on | out_name: Name for output dir instead of the current data
- hf | hand_feat: If you want to start from hand curated feature files meaning starting already with pos and neg input file
- j | n_jobs: Number of jobs used for graph feature computation


#### Output 
At the end of the run the location of the result table is given.
The final result table will have all coulms of the input table and the additional a prediction column, where you find the predicted class of each RRI (0 or 1).

Thought the program seveal output files are generated in the following structure:

    ├── date_Cherri\_evaluation\_mode
    |   ├── evaluate\_RRIs.table
    |   ├── positive\_instance
    |       ├── test\_eval\_context\_{context}pos.csv
    |       ├── date_occ_out
    |           ├── occupied\_regions.obj
    |           ├── rri_occupied\_regions\_overlapTH\_0.3\_scoreTH\_1.csv
    |   ├── feature\_files
    |       ├── feature\_filtered\_test\_eval\_context\_150\_pos.csv
    |       ├── training\_data\_test\_eval\_context\_150.npz
    |   ├── evaluation
    |       ├── evaluation\_results\_test\_eval.csv


### Second mode **train** : build new Cherri model
To train mode of Cherri enables you to train your own model. 
Using Chira RRI output data 'ChiRA interaction summary' table for "build model" mode will generate a prediction classifier.
If you you want to prepare the 'ChiRA interaction summary' table please follow the [tutorial](https://training.galaxyproject.org/training-material//topics/transcriptomics/tutorials/rna-interactome/tutorial.html).


#### Input format of the RRIs to evaluate
Please select the Chira RRI output files as input for cherri train. You should specify the path to the folders containing all replicats. If you want to filer for interactions, where Chira also alrady found a hyprid within its analysis you should enalble the hybritsation with the Chira workflow and set for Cherri the filter_hybrid pararmter as 'on'. Some of the possible interaction maid be missed this way, but it will also filter out potential false postive interactions. 


#### Input Parameter
##### required
- i1 | RRI_path: Path to folder storing the ChiRA interaction summary files
- o | out_path: Path to folder all output folder of each step of the data preparation
- r | list\_of\_replicates: List ChiRA interaction summary files names of all replicates
- l | chrom\_len\_file: Tabular file containing data in two columns one row for each chromosome: 'chrom name' \t 'chrom lenght'
- g | genome: Path to 2bit genome file
##### optional
- c | context: How much context should be added at left an right of the sequence
- n | experiment_name: E.g. name of the data source of RRIs. Will be used for the file names
- p | param\_file: IntaRNA parameter file
- st | use_structure: Set 'off' if you want to disable graph-kernel features default: 'on' (when set to 'on' the feature optimization will be performed directly and the data will be stored in feature_files and no model/feature folder will be set up)
- i2 | RBP_path: Path to the genomics RBP cross link positions location (in bed format)
- t | run_time: Time used for the optimization in seconds default: 12h
- me | memoryPerThread: Memory in MB which each thread can use (total ram/threads)
- j | n_jobs: Number of jobs for the optimization
- mi | mixed: Use mixed model to combine different dataset into a combined model. You need to add all precomputed models, who's data you would like to combine in one single folder. Rename the first level folder to the experiment name and speciy the path to this folder as i1 parameter and the experiment names as replicates (-r). The Cherri pipeline will than first concatenate all positive and negative data set file of the given experiments and stars the Cherri call form the feature generation step. Defualt. 'off' (advanced option) 
- fh |filter_hybrid: Filter the data for hybrids already detected by chira (set to 'on' to filter default:'off')



#### Output 
At the end of the run the location of the trained model is given.

Thought the program seveal output files are generated in the following structure:

    ├── date_Cherri_model_build
    |   ├── date_occ_out
    |       ├── occupied_regions.obj
    |       ├── rri_occupied_regions_overlapTH_0.3_scoreTH_1.cvs
    |   ├── read_pos_neg_data
    |       ├── test_train_context_50_pos_occ_neg.csv
    |       ├── test_train_context_50_pos_occ_pos.csv
    |   ├── feature_files
    |       ├── feature_filtered_test_eval_context_150_pos.csv
    |       ├── feature_filtered_test_eval_context_150_neg.csv
    |       ├── training_data_test_eval_context_150.npz (filtered features if use_structure==on)
    |   ├── model
    |       ├── features
    |           ├── test_train_context_50.npz (only present when use_structure==off)
    |       ├── optimized
    |           ├── test_train_context_50.model
    |           ├── test_train_context_50.cvs




### Get to know Cherri with a test run 
Using test data stored in this reposetory under 'test_data' you can see how Cherri is running. First please install Cherri via conda. 

#### Test Cherri train model mode
```
cherri eval -i1 /vol/scratch/Cherri/test_data/evaluate/test_evalueat_rris.cvs -g /vol/scratch/data_storage/data/genomes/hg38_UCSC_20210318.2bit -l /vol/scratch/data_storage/data/genomes/hg38_Info.tab -o /vol/scratch/data_storage/ -n test_eval -c 150 -st on
```

#### Test Cherri train evaluation mode
```
cherri train -i1 ./Cherri/test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g /vol/scratch/data_storage/data/genomes/hg38_UCSC_20210318.2bit -l /vol/scratch/data_storage/data/genomes/hg38_Info.tab -o /vol/scratch/data_storage/ -n paris_human_test -c 50 -st on -t 10000 -me 8000 -j 7
```

### Backround of scripts which are wrapped in Cherri
Cherri is build as a modular tool calling individual scrips doing specific tasks. If you would like to only perform one step of the Cherri piplein you can do this by calling the individual scrips. A short description of this scripts is given in the following.


### find_trusted_RRI.py
Here we search for trusted RRIs, so RRIs which can be found in all replicates. In a first filter step only uniquely mapped RRIs are taken. Than RRI sequence partners in all replicas are found, using a overlap threshold. Output are the Chira input tables, now containing only the trusted RRIs. Only one of the sequence RRI pairs is added to the output. 

#### example call
```
python find_trusted_RRI.py -i /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/data/data_Chira/training/Paris -r test_rep1.tabular test_rep2.tabular test_rep3.tabular -o 0.6 -d /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/ -n test_paris
```

#### Input Parameter
- input_path: path to folder storing all input data
- list_of_replicates: list of replicates
- overlap_th: path output repository
- experiment_name: name of the data source of positive trusted RRIs

#### Output 
- trusted RRIs in tabular format
- pickled list of avg overlap percent
- pickled list of avg overlap length


### plot_avg_overlaps.py 
The scirpt takes the lists of avg overlap percent and overlap length form the find_trusted_RRI.py script. The Idea is that for different overlap thresholds will be plotted into a box plot. At the moment for overlap thresholds ['0.3', '0.4', '0.5', '0.6', '0.7']. Therefore, the find_trusted_RRI.py should be called with this overlap thresholds.

#### example call
```
python plot_avg_overlaps.py -i /vol/scratch/data/plots/overlap_distibution/ -n rri_overlap_plots
```

#### Input Parameter
- input_path: path to folder where input files
- experiment_name: name of the data source of positive trusted RRIs


#### Output 
- two box plots



### find_occupied_regions.py
Given the RRI information tables form Chira and RNA-Protein binding positions a Interlab object is build. The occuped information can be used to mask parts of the genome and therefore enable to select negative interaction regions. This regions are not part of interaction in nature. 

#### example call
```
python find_occupied_regions.py -i1 /vol/scratch/data/RRIs/Paris/ -r test_rep1.tabular test_rep2.tabular test_rep3.tabular -o /vol/scratch/data/RRIs/
```

#### Input Parameter
- RRI_path: path to folder storing all RRI data (tabular)
- rbp_path: path to RBP side data file (bed format)
- list_of_replicates: list having filenames of all replicates
- out_path: path to folder storing outputfiles
- overlap_th: overlap threshold
- score_th: score threshold

#### Output 
- outputs a pickeled Interlab object and prints the path to the file
- since it calls the trusted RRIs script the output will also be a file for trusted RRIs once with a lower score threshold to build the occupied regions from and a score filtered file, which can be used to build the build positive and negative prediction instances.


### generate_pos_neg_with_context.py
Given trusted RRI and the occupied regions a given context is appended. Than positive sequences are computed by calling IntaRNA specifying the the seed region with the trusted RRI regions. The negative interactions are computed by calling IntaRNA given regions, which should not be in the interaction 'occupied regions'. 



#### example call
```
python generate_pos_neg_with_context.py -i1 /vol/scratch/data/trusted_RRIs/test_paris_overlap_0.6.cvs -i2 /vol/scratch/data/RRIs/20210809_occ_out/occupied_regions.obj -d /vol/scratch/data/pos_neg_data_context/ -g /vol/scratch/data/genomes/hg38_UCSC_20210318.2bit -n test_paris_HEK293T  -c 10
```




#### Input Parameter
- input_rris: path to file storing all positve trusted RRIs
- input_occupyed: path to file storing occupied sides (InterLab obj)
- output_path: path output repository
- experiment_name: name of the data source of positive trusted RRIs
- genome_file: path to 2bit genome file
- context: how much context should be added at left an right of the sequence
- pos_occ/ no_pos_occ: if positive instances should have the occupied positions blocked with the IntaRNA call
- block_ends: number of blocked nucleotides at the start and end of the sequence
- no_sub_opt: # of interactions IntraRNA will give is possible
- chrom_len_file: tabular file containing Chromosome name \t Chromosome length for each chromosome
- param_file: IntaRNA parameter file


#### IntaRNA params for call

To generate the current features we need:

| param  | value  | description | 
|---|---|---|
| outMode |  C  | Output style (C=tabluar format) |
| seedBP |  5  | the number of base pairs within the seed |
| seedMinPu |  0  | the minimal unpaired probability of each seed region in query and target |
| accW |  150  | sliding window length (0=global folding) |
| acc |  N/C  |  To globally turn off accessibility consideration: turn off/on |
| outMaxE |  -5  | maximal energy for any interaction reported |
| outOverlap |  B  | overlapping of interaction sites of suboptimal allowed (B:both) |
| outNumber |  5  | generate up to N interactions for each query-target pair |
| seedT/QRange |  positive interaction | genomic positions of the trusted RRI |
| q/tAccConstr |  negative interaction  | genomic positions of the occupied regions |
| intLenMax |  50 | restrict the overall length an interaction |
| temperature |  37  | experimental temperature |
| intLoopMax |  3  | number of unpaired bases between inter molecular base pairs |



#### Interaction prediction oputput

To generate the current features we need:

| description  | coloum name  | 
|---|---|
| start position of RRI target  | 'start1'  |
| end position of RRI target   |  'end1' | 
|  start position of RRI query  |  'start2' |
| end position of RRI query   |  'end2' | 
| dotbrackad of interaction  | hybridDP  |
| RRI sequence target&query  |  subseqDP | 
| MFE  |  E | 
| seeds starts listed for target |  seedStart1 | 
| interacting target sequence  |  target | 
| interacting query sequence | query |


#### Output 
- Positive and negative datasets stored in tabular format. 


### plot_tRRIs.py
For the positive and negative datasets overview plots are generated. 

#### example call
```
python plot_tRRIs.py -i1 test_paris_HEK293T_context_method_together_shuffling_method_3_pos_RRI_dataset.csv -i2 test_paris_HEK293T_context_method_together_shuffling_method_3_neg_RRI_dataset.csv -i3 /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/ -o /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/
```

#### Input Parameter
- input_pos: file name of positive dataset in tabular format (, separted)
- input_neg: file name of negative dataset in tabular format (, separted)
- save_path: directory where the plots will be stored in a plot folder

#### Output 
- RNA histogram plot
- Energy distribution plot


### get_features.py
Here additional sequence features are computed and the output is filtered for the final or a given feature set. The features are stored in a tabular format.



#### example call
```
python get_features.py -i ../output/paris_HEK293T_06_context_method_together_shuffling_method_2_pos_RRI_dataset.csv -f E no_bps GC_content mfe_normby_GC -o /home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/input_features/2
```

#### Input Parameter
- input: path to input file
- feature_set_list: list of all features that should be summarized in the output
- output_file: file path where the output table should be stored

#### Output 
- tabular file having all features given via the feature_set_list


### feature selection and optimization

If you would only like to run the feature selection or optimization please check out the [biofilm](https://github.com/smautner/biofilm). 


