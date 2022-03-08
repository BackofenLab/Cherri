# :cherries: Cherri (Computational Help Evaluating RNA-RNA interactions) :cherries:

## Progamm Idea:
Cherri is a tool that distinguishes between biological relevant or real RRIs and RRIs which most likely would not occurs in nature. Giving the tool a RRI prediction it will evaluate weather a prediction is likely to be a real prediction. 
Software we developed models base on the human and mouse RNA-RNA interactome. 
Further more the option of training a novel model basted on RNA-RNA interactome data is given.


## Instalation

Cherri developed in Linux and tested on Ubuntu (18.04 LTS). For the installation on your conda is requierd. 


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the Python 3.8 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be evocable on the command line via (possibly in a more advanced version):

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
Or install it into your favorite existing enviorment
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

conda install -c conda-forge scikit-learn
conda install -c conda-forge networkx
conda install -c bioconda ucsc-twobittofa
conda install -c bioconda interlap
conda install -c bioconda IntaRNA
conda install -c bioconda numpy
conda install -c bioconda pandas
conda install -c smautner eden-kernel
conda install -c smautner biofilm=0.0.78
```
Or directly:

```
conda create -n cherri -c conda-forge -c bioconda -c smautner scikit-learn networkx ucsc-twobittofa interlap pandas IntaRNA eden-kernel biofilm=0.0.78

```

Cherri was tested with the following versions: XX. 
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



## Usage
You can use Cherri in two modes. The first mode to predict weather your RRI prediction biological relevant. You can used the current model for human or mouse data. If you would like build a model based on novel RRI interactome data you can use the second mode.

### First mode **eval** :  evaluation of RRIs
Based on a tabular file containing chromosomal position data of the RRIs, it classify if the interaction region is likely to be a relevant one.

Cherri has a model traind on human data integrated. Cherris model can be used to predict RRIs of different organisms. However, if there exists a RNA-RNA-interactome dataset for your preferred organism we recommend to train your own organism specific model using cherris train mode. This model can be than use in the eval mode for the classification of your predicted RRIs.

#### Input format of the RRIs to evaluate
The instances which should be tested should be in a tabular format. You should specify a header lines, with the chromosome interaction start interaction end and strand of the two interacting partners. 
```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```
With a specific setting only 'positive' instances are computed. If no occupied regions are not given by using none this can be specified. If the occupied regions for human or the mouse model used for training should be used specify mouse or human. 


#### Input Parameter
- i1 | RRIs_table: table containing all RRIs that should be evaluated in the correct format
- g | genome_file: path to 2bit genome file
- o | out_path: path to folder all output folder of each step of the data preparation
- l | chrom_len_file: tabular file containing data in two columes one row for each chromosome: 'chrom name' \t 'chrom lenght' 
- i2 | occupyed_regions: path to occupied regions object file. give human, mouse, non or your own file path.
- c | context: how much context should be added at left an right of the sequence
- n | experiment_name: name of the data source of RRIs
- p | param_file: IntaRNA parameter file
- m | model_file: set if a model created by train module should be used
- mp | model_params: set path to feature file of new model if model_file is changed
- st | use_structure: set 'off' if you want to disable structure, default 'on'


#### Output 
At the end of the run the location of the result table is given. 
The final result table will be like the input table and the additional prediction column, where you find the class of the RRI (0 or 1).

Thought the program seveal output files are generated in the following structure:

    ├── date_Cherri_evaluation_mode
    |   ├── evaluate_RRIs.table
    |   ├── positive_instance
    |       ├── test_eval_context_{context}pos.csv
    |       ├── date_occ_out
    |           ├── occupied_regions.obj
    |           ├── rri_occupied_regions_overlapTH_0.3_scoreTH_1.cvs
    |   ├── feature_files
    |       ├── feature_filtered_test_eval_context_150_pos.csv
    |       ├── training_data_test_eval_context_150.npz
    |   ├── evaluation
    |       ├── evaluation_results


### Second mode **train** : build new Cherri model
Using Chira RRI output data 'ChiRA interaction summary' table for "build model" mode will generate a prediction classifier. 


#### Input format of the RRIs to evaluate
Please select the Chira RRI output files as input for cherri train. 


#### Input Parameter
- i1 | RRI_path: path to folder storing the ChiRA interaction summary files
- o | out_path: path to folder all output folder of each step of the data preparation
- r | list_of_replicats: list ChiRA interaction summary files names of all replicates
- l | chrom_len_file: tabular file containing data in two columes one row for each chromosome: 'chrom name' \t 'chrom lenght'
- g | genome_file: path to 2bit genome file
- c | context: how much context should be added at left an right of the sequence
- n | experiment_name: name of the data source of RRIs
- p | param_file: IntaRNA parameter file
- st | use_structure: set off if you want to disable structure 
- i2 | path to RBP binding location file in bed format


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
    |       ├── training_data_test_eval_context_150.npz
    |   ├── model
    |       ├── features
    |           ├── test_train_context_50.npz
    |       ├── optimized
    |           ├── test_train_context_50.model
    |           ├── test_train_context_50.cvs





### Backround of scripts which are wrapped in Cherri
The following scripts are needed to build the trainings data to build a model.


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
The scirpt takes the lists of avg overlap percent and overlap lenght form the find_trusted_RRI.py script. Thie Idea ist that for differen overlap thresholds will be plotted into a box plot. At the moment for overlap thresholds ['0.3', '0.4', '0.5', '0.6', '0.7']. Therefore, the find_trusted_RRI.py should be called with this overlap thresholds.

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
- outputs a pickeled Interlab object and prints the path to the file
- since it calls the trusted RRIs script the output will also be a file for trusted RRIs once with a lower score threshold to build the occupied regions from and a score filtered file, which can be used to build the build positive and negative prediction instances.


### generate_pos_neg_with_context.pl
Given trusted RRI and the occupied regions a given context is appended. Than positive sequences are computed by calling IntaRNA specifiying the the seed region with the trusted RRI regions. The negative interactions are computed by calling IntaRNA given regions, which should not be in the interaction 'occupied regions'. 



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
- pos_occ/ no_pos_occ: if postived instances should have the occupied positions blocked with the IntaRNA call
- block_ends: number of blocked nucleotides at the start and end of the sequence
- no_sub_opt: # of ineractions IntraRNA will give is possible
- chrom_len_file: tabular file containing chrom name \t chrom lenght for each chromosome
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
| outOverlap |  B  | overlapping of interaction sites of Suboptimal allowed (B:both) |
| outNumber |  5  | generate up to N interactions for each query-target pair |
| seedT/QRange |  positive interaction |  |
| q/tAccConstr |  negative interaction  |  |
| intLenMax |  50 | restrict the overall length an interaction |
| temperature |  37  | experimental temperature |
| intLoopMax |  3  | number of unpaired bases between intermolecular base pairs |



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
Here additional sequence features are computed and the output is filtered for the final or a given feature set. The features are stored in a tabular format.


#### Our List of features
- E : Minimum free energy (overall interaction energy)
- E_hybrid : energy of hybridization only = E - ED1 - ED2
- maxED: maximal energy of ED1 and ED2
- no_bps : number of base pairs within the interaction
- GC_content: GC content of the interaction side
- max_inter_len: the maximum interaction side calculated by the length of the target sequence and query sequence length
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


