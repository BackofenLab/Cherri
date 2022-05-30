# :cherries: CheRRI (Computational Help Evaluating RNA-RNA interactions) :cherries:

Cherri detects functional RNA-RNA interactions (RRI) sites, by evaluating if a interaction site most likely occurs in nature.
It will help to filter false positive interaction sites, which can be generated either experimentally or by a RRI prediction algorithms. 

## Progamm Idea:
CheRRI can be run in two modes, the model generation or **train** mode, or the RRI evaluation or **eval** mode. Here is a illustration of CheRRI's workflow:

<img src="./plots/Cherri_workflow_resuctured2.svg " alt="Cherri_workflow" width="550"/>


For the evaluation of given RRI sites a model must be specified in Cherri's **eval** mode. Here pre-trained models can be applied or the user trains a model it self using CheRRI's **train** mode.
To train a novel model RNA-RNA interactome specifying all RRIs sites within a organism should be provided. CheRRI makes use of replicate data by checking if a RRI site can be found in all replicates within a overlap threshold. This is how CheRRI builds the set of trusted RRIs. In the evaluation mode the interaction positions are reformatted in the same way as the trusted RRIs. Than in both cases CheRRI's core method can be run to generate a feature set which will be used to select a model in the **train** and in the **eval** mode used for the evaluation of the biological relevance for the submitted RRI sites.




## Instalation

CheRRI is developed in Linux and tested on Ubuntu (18.04 LTS). For the installation on your Conda is required. 


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the Python 3.8 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be accessed on the command line via (possibly in a more advanced version):

```
$ conda --version
conda 4.10.3
```

### Install CheRRI Conda package **not working jet**

CheRRI is available as a Conda package [here](https://anaconda.org/bioconda/).

We recommend to create a new Conda environment inside which we will then install CheRRI:

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

You additionally need to fix the python hash seed:

```
conda env config vars set PYTHONHASHSEED=31337
```
Please reactivate the environment to activate the PYTHONHASHSEED environment variable.
Now CheRRI should be available and usable inside the environment:


```
cherri -h
```

### Manual installation

To manually install CheRRI, you first create a Conda environment. Inside the environment, you need to install the following dependencies:

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

Finally, to install the tool itself, you simply clone the repository and execute the installation script inside the folder:

```
git clone https://github.com/BackofenLab/Cherri.git
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

Now we can run CheRRI from any given folder (just remember to re-activate the environment once you open a new shell):

```
cherri -h
```

CheRRI was tested with the following versions: 0.1. 

## Usage
You can use CheRRI in two modes. The eval mode predicts weather the site of your RRI prediction is biological relevant. Pre-computed models exist for human or mouse data and can be downloaded. If you would like to build a model based on novel RRI interactome data you can use the **train** mode.
If you cloned CheRRI's repository from GitHub you can apply CheRRI example calls. All example the following example calls should be working if you are one level before the CheRRI folder and have the conda environment activated.

### Mode **eval**:  evaluation of RRIs
Based on a tabular file containing chromosomal position data of the RRIs, CheRRI classify if the interaction region is likely to be a biological relevant one.

For the **eval** mode a model and the filtered features set have to be specified.
CheRRI's has pre-trained models for human and mouse. 
The pre-trained models can be downloaded from [zenodo](https://zenodo.org/record/6533932#.Ynu--FxBwUE)  DOI 10.5281/zenodo.6533931.
If there exists a RNA-RNA-interactome dataset for your preferred organism we recommend to train your own organism specific model using CheRRI's **train** mode. After the training this model can be than use, here in the **eval** mode, for the classification of your predicted RRIs positions.


#### Input format of the RRI evaluation table (RRIs_table)
The to evaluate instances should in a tabular format. You should specify a header lines, with the chromosome interaction start interaction end and strand of the two interacting partners:
```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```
With a specific setting only 'positive' instances are computed. If no additional occupied regions are specified only the once of the given input interactions are used. However, we recommend to specify the occupied regions of the trained model, to be consisted within the feature generation.

#### Example call CheRRI **eval** mode
For the test call please download the [Cherri_models_data](https://zenodo.org/record/6533932#.Ynu--FxBwUE).
```
cherri eval -i1 /vol/scratch/Cherri/test_data/evaluate/test_evalueat_rris.cvs -g human -l human -o ./ -n test_eval -c 150 -st on -m ./Cherri_models_data/Model_with_graph_features/PARIS_human/model/full_PARIS_human_context_150.model -mp ./Cherri_models_data/Model_with_graph_features/PARIS_human/feature_files/training_data_PARIS_human_context_150.npz
```



#### Input Parameter
##### required:
| ID | name | description |
|---|---|-----|
| `-i1` |`--RRIs_table` | Table containing all RRIs that should be evaluated in the correct format|
| `-g` | `--genome_file`| Path to 2bit genome file, or used the build in download if you want the human or moues genomes |
| `-o` | `--out_path`| Path to folder all output folder of each step of the data preparation |
| `-l` | `--chrom_len_file` | Tabular file containing data in two columns one row for each chromosome: 'chrom name' \t 'chrom length'. You can specify 'human' or 'mouse' directly. |
| `-m` | `--model_file` | Set path to the model which should be used for evaluation |
| `-mp` | `--model_params` | Set path to feature file of the given model |
##### optional:
| ID | name | description |
|---|---|-----|
| `-i2` | `--occupyed_regions` | Path to occupied regions object file |
| `-c` | `--context` | How much context should be added at left an right of the sequence |
| `-n` | `--experiment_name` | Name of the data source of RRIs |
| `-p` | `--param_file` | IntaRNA parameter file. default: file in rrieval/IntaRNA_param |
| `-st` | `--use_structure` | Set 'off' if you want to disable structure, default 'on' |
| `-on` | `--out_name` | Name for output dir instead of the current data |
| `-hf` | `--hand_feat` | If you want to start from hand curated feature files meaning starting already with pos and neg input file |
| `-j` | `--n_jobs` | Number of jobs used for graph feature computation |


#### Output 
At the end of the run the location of the result table is given.
The final result table will have all columns of the input table and the additional a prediction column, where you find the predicted class of each RRI (0 or 1).

Thought the program seveal output files are generated in the following structure:

    ├── date_Cherri_evaluation_mode
    |   ├── evaluate_RRIs.table
    |   ├── positive_instance
    |       ├── test_eval_context_{context}pos.csv
    |       ├── date_occ_out
    |           ├── occupied_regions.obj
    |           ├── rri_occupied_regions_overlapTH_0.3_scoreTH_1.csv
    |   ├── feature_files
    |       ├── feature_filtered_test_eval_context_150_pos.csv
    |       ├── training_data_test_eval_context_150.npz
    |   ├── evaluation
    |       ├── evaluation_results_test_eval.csv


### Mode train: built new CheRRI model
Within CheRRI's train mode you to train your own model. 
Using ChiRA RRI output data 'ChiRA interaction summary' table for "build model" mode will generate a prediction classifier.
If you want to prepare the 'ChiRA interaction summary' table please follow the [tutorial](https://training.galaxyproject.org/training-material//topics/transcriptomics/tutorials/rna-interactome/tutorial.html).


#### Retive RNA-RNA interactome files
Please select the ChiRA RRI output files as input for CheRRI train mode as replicates. You should specify the path to the folders containing all replicates. If you want to filter for interactions, where ChiRA already found a hybrid within its analysis you should enable the hybridization with the ChiRA workflow and set within CheRRI's input parameters the filter_hybrid parameter to 'on'. Some of the possible interaction maid be missed this way, but it will also filter out potential false positive interactions. 

#### Example call CheRRI train model mode
```
cherri train -i1 ./Cherri/test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7
```



#### Input Parameter
##### required
| ID | name | desctiption |
|---|---|-----|
| `-i1` | `--RRI_path`| Path to folder storing the ChiRA interaction summary files|
| `-o` | `--out_path`| Path to folder all output folder of each step of the data preparation |
| `-r` | `--list_of_replicates`| List ChiRA interaction summary files names of all replicates |
| `-l` | `--chrom_len_file`| Tabular file containing data in two columns one row for each chromosome: 'chrom name' \t 'chrom lenght' |
| `-g` | `--genome`| Path to 2bit genome file
##### optional
| ID | name | desctiption |
|---|---|-----|
| `-c` | `--context`| How much context should be added at left an right of the sequence |
| `-n` | `--experiment_name`| E.g. name of the data source of RRIs. Will be used for the file names|
| `-p` | `--param_file`| IntaRNA parameter file |
| `-st` | `--use_structure`| Set 'off' if you want to disable graph-kernel features default: 'on' (when set to 'on' the feature optimization will be performed directly and the data will be stored in feature_files and no model/feature folder will be set up)|
| `-i2` | `--RBP_path`| Path to the genomics RBP cross link positions location (in bed format) |
| `-t` | `--run_time`| Time used for the optimization in seconds default: 12h|
| `-me` | `--memoryPerThread`| Memory in MB which each thread can use (total ram/threads)|
| `-j` | `--n_jobs`| Number of jobs for the optimization|
| `-mi` | `--mixed`| Use mixed model to combine different dataset into a combined model. Default. 'off' (advanced option)| 
| `-fh` |`--filter_hybrid`| Filter the data for hybrids already detected by ChiRA (set to 'on' to filter default:'off') |

For the mixed mode parameter you need to add all pre-computed models, who's data you would like to combine in one single folder. Rename the first level folder to the experiment name and specify the path to this folder as i1 parameter and the experiment names as replicates (-r). The CheRRI pipeline will than first concatenate all positive and negative data set file of the given experiments and stars the CheRRI call form the feature generation step. 


#### Output 
At the end of the run the location of the trained model is given.

Thought the program several output files are generated in the following structure:

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




## CheRRI's core method scripts
CheRRI is build as a modular tool calling individual scrips doing the tasks of the CheRRI's core method. If you would like to only perform one step of the CheRRI pipeline you can do this by calling the individual scrips. A short description of this scripts is given in the following.


### RRI detection: find_trusted_RRI.py
Here we search for trusted RRIs, so RRIs which can be found in all replicates. In a first filter step only uniquely mapped RRIs are taken. Than RRI sequence partners in all replicas are found, using a overlap threshold. Output are the ChiRA input tables, now containing only the trusted RRIs. Only one of the sequence RRI pairs is added to the output. 

#### Input Parameter
| ID | name | description |
|---|---|-----|
| `-i` | `--input_path` |  Path to folder storing input data (replicates) |
|`-r`| `--list_of_replicats` | Filenames list of all replicates |
| `-o` | `--overlap_th` | Overlap threshold to find trusted RRIs |
| `-d` | `--output_path` | Path where output folder should be stored |
|`-n` | `--experiment_name` | Name of the data soruce of positive trusted RRIs |
| `-s` | `--score_th` | Threshold for EM score from ChiRA |
| `-fh` | `--filter_hybrind` | Filter the data for hyprids alrady detected by ChiRA

#### Output 
The filtered set of trusted RRI sites in tabular format. 

### Compute occupied regions: find_occupied_regions.py
Given the RRI information tables form ChiRA and RNA-Protein binding positions a InterLab object is build. The occupied information can be used to mask parts of the genome and therefore enable to select negative interaction regions. This regions are not part of interaction in nature. 


#### Input Parameter
| ID | name | description |
|---|---|-----|
|`-i1` | `--RRI_path` | Path to folder storing all RRI data (tabel) |
| `-i2` | `--rbp_path` | Path to RBP side data file (BED format) |
| `-r` | `--list_of_replicats` | Filenames list of all replicates |
| `-o` | `--out_path` | Path where output folder should be stored |
| `-t` | `--overlap_th` | Overlap threshold |
| `-s` | `--score_th` | Score threshold | 
| `-e` | `--external_object` |External RRI  overlapping object (InterLap dict)|
| `-fh` | `--filter_hybrind` | Filter the data for hybrids already detected by ChiRA|
| `-mo` | `--mode` | Function call within which CheRRI mode [train/eval]|

#### Output 
A python pickled file object storing occupied regions in a InterLap directory. 


### Interaction predictions: generate_pos_neg_with_context.py
Given trusted RRI and the occupied regions a given context is appended. Than positive sequences are computed by calling IntaRNA specifying the the seed region with the trusted RRI regions. The negative interactions are computed by calling IntaRNA given regions, which should not be in the interaction 'occupied regions'. 


#### Input Parameter

| ID | name | description |
|---|---|-----|
| `-i1` | `--input_rris` |Path to file storing all trusted RRIs|
|  `-i2` | `--input_occupyed` | Path to occupying regions file |
| `-d` | `--output_path` | Path where output folder should be stored |
| `-n` | `--experiment_name` |Name of the data source of positive trusted RRIs|
| `-g` | `--genome_file` | Path to 2bit genome file |
| `-c` | `--context` |How much context should be added up- and downstream |
|   | `--pos_occ`  | Default occupied regions are set |
|   | `--no_pos_occ`  | Set if no occupied regions should be set |
| `-b` | `--block_ends` |# nts blocked at the ends of the sequence|
| `-s` | `--no_sub_opt` | # of interactions IntraRNA will give is possible |
| `-l` | `--chrom_len_file` |Tabular file containing chrom name \t chrom length for each chromosome |
| `-p` | `--param_file` | IntaRNA parameter file |
| `-m` | `--mode` | Which CheRRI mode is running|



#### IntaRNA params for call

To generate the current features IntaRNA parameters are set link this:

| parameters  | value  | description | 
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

IntaRNA parameters can be changed by specifying a IntaRNA parameter file.



#### Output 
Positive and negative datasets stored in tabular format. 


### Feature extraction: get_features.py
Here additional sequence features are computed and the output is filtered for the final or a given feature set. The features are stored in a tabular format.
Note that the graph-kernel features are computed directly in CheRRI's main functions and do not have a sparat file. 

#### Input Parameter
| ID | name | description |
|---|---|-----|
| `-i` | `--input` | Path to input file |
| `-f` | `--feature_set_list` | Set of features the script will output |
| `-o` | `--output_file` | Output file path inclusive of the file name |


#### Prediction output used to build the features

The following columns are used to build the set of features:

| description  | column name  | 
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
- tabular file having all features given via the feature_set_list


### Feature selection and optimization

If you would only like to run the feature selection or optimization please check out the [biofilm](https://github.com/smautner/biofilm). 



