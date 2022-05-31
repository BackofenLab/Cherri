# :cherries: CheRRI (Computational Help Evaluating RNA-RNA interactions) :cherries:

CheRRI detects functional RNA-RNA interaction (RRI) sites, by evaluating if an interaction site most likely occurs in nature.
It helps to filter interaction sites generated either experimentally or by an RRI prediction algorithm, by removing false positive interactions.



## CheRRI's workflow

CheRRI can be run in two modes, the model generation or **train** mode, or the RRI evaluation or **eval** mode. Here is an illustration of CheRRI's workflow:

<img src="./plots/Cherri_workflow_resuctured2.svg " alt="Cherri_workflow" width="550"/>


For the evaluation of a given set of RRI sites, a model must be specified in CheRRI's **eval** mode. Here pre-trained models can be applied or the user trains a model using CheRRI's **train** mode.
To train a novel model, an RNA-RNA interactome dataset specifying all RRI sites should be provided. CheRRI makes use of replicate data by checking if an RRI site can be found in all replicates within an overlap threshold. This is how CheRRI builds the set of trusted RRIs. In **eval** mode, the interaction positions are reformatted in the same way as the trusted RRIs. In both modes, CheRRI uses the same core method to generate a feature set which will be used to select a model in **train** mode, and in **eval** mode for the evaluation of the biological relevance of the submitted RRI sites.




## Installation

CheRRI was developed in Linux and tested on Ubuntu (18.04 LTS). Conda is required to install CheRRI.


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the newest Python 3 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be accessed on the command line via (note that your version can be different):

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
Alternatively, you can install it into your existing environment via:

```
conda install -c bioconda cherri
```

If you are experiencing problems while running `conda install -c bioconda cherri` (e.g. complaints about conflicting dependencies), the following commands might help:

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

This tells Conda to explicitly look for packages in the specified channels, stored in the `.condarc` [conda configuration file](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html).

You additionally need to fix the python hash seed:

```
conda env config vars set PYTHONHASHSEED=31337
```
Please reactivate the environment to activate the PYTHONHASHSEED environment variable.

```
conda deactivate
conda acivate cherri
```

Now CheRRI should be available and usable inside the environment:


```
cherri -h
```

### Manual installation

To manually install CheRRI, first create a Conda environment:

```
conda create -n cherri python=3.8 -c conda-forge -c bioconda
conda activate cherri
```
Inside the environment, you need to install the following dependencies:


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

Or create the environment with all dependencies at once:

```
conda create -n cherri -c conda-forge -c bioconda -c smautner scikit-learn networkx ucsc-twobittofa interlap pandas IntaRNA python-wget eden-kernel biofilm=0.0.88 python-wget
conda activate cherri
```

Finally, to install the tool itself, simply clone the repository and execute the installation script inside the cloned folder:

```
git clone https://github.com/BackofenLab/Cherri.git
cd Cherri
python -m pip install . --ignore-installed --no-deps -vv 
```

You additionally need to set a fixed python hash seed within the conda environment:

```
conda env config vars set PYTHONHASHSEED=31337
```
Or set it just for your current session:

```
export PYTHONHASHSEED=31337
```
After setting the environment variable, reactivate your environment:
```
conda deactivate
conda acivate cherri
```

Now you can run CheRRI from any given folder:

```
cherri -h
```



## Usage

You can use CheRRI in two modes. The **eval** mode predicts whether an RRI site of interest is biologically relevant. Pre-computed human and mouse models exist and can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.6533932). If you would like to build a model based on novel RRI interactome data, you can use the **train** mode.

### Evaluation of RRIs

Based on a tabular file containing chromosomal position data of the RRIs, CheRRI classifies if the interaction region is likely to be a biologically relevant one.

For the **eval** mode, a model and the filtered feature set have to be specified.
CheRRI has pre-trained models for human and mouse, which can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.6533932) (see `content.txt` inside the zip folder for details).
If there exists an RNA-RNA-interactome dataset for your preferred organism, we recommend to train your own organism-specific model using CheRRI's **train** mode. After training, the model can be than used in **eval** mode for the classification of your predicted RRI positions.


#### RRI input format in evaluation mode

The RRI instances to be evaluated need to be given in tabular format (parameter `--RRIs_table`). The table needs the following header line:

```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```

Following the header line, each subsequent line represents an RRI, with chromosome ID (format: 1,2,3 ...), interaction start, interaction end, and strand ("+" or "-") of the two interacting partners. For example, you might want to evaluate the following three RRI sites:

```
19,18307518,18307539,-,14,90454500,90454521,+
X,109054541,109054590,+,9,89178539,89178562,-
10,123136102,123136122,+,5,1245880,1245902,+
```

If no additional occupied regions are specified (`--occupied_regions`), only the ones of the given input interactions (`--RRIs_table`) are used. However, to be consistent with the feature generation of the trained model, we recommend that the user also specifies the occupied regions used to train the model provided via `--model_file`. For example, for the PARIS_human model without graph features, the occupied regions data object is located at 
`Model_without_graph_features/PARIS_human/occupied_regions/occupied_regions.obj` (inside the mentioned zip file from [Zenodo](https://doi.org/10.5281/zenodo.6533932)).





#### Example call CheRRI **eval** mode

For the test call please download the [Cherri_models_data](https://doi.org/10.5281/zenodo.6533932) zip folder. The PARIS_human model is needed to execute the call. Be sure to provide the correct location for the model and its feature set (`-m`, `-mp`). For example, assuming the data (zip folder extracted to folder `Cherri_models_data`) is stored inside the Cherri folder:

```
cherri eval -i1 test_data/evaluate/test_evaluate_rris.cvs -g human -l human -o ./ -n test_eval -c 150 -st on -m Cherri_models_data/Model_with_graph_features/PARIS_human/model/full_PARIS_human_context_150.model -mp Cherri_models_data/Model_with_graph_features/PARIS_human/feature_files/training_data_PARIS_human_context_150.npz
```



#### Input Parameters

Input parameters for Cherri's **eval** mode (`cherri eval`):

##### required:
| ID | name | description |
|---|---|-----|
| `-i1` |`--RRIs_table` | Table containing all RRIs that should be evaluated in the correct format|
| `-g` | `--genome_file`| Path to 2bit genome file, or used the build in download if you want the human or mouse genome |
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data and feature preparation as well as the evaluated instances |
| `-l` | `--chrom_len_file` | Tabular file containing data in two columns format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-m` | `--model_file` | Set path to the model which should be used for evaluation |
| `-mp` | `--model_params` | Set path to the feature file of the given model |
##### optional:
| ID | name | description |
|---|---|-----|
| `-i2` | `--occupied_regions` | Path to occupied regions python object file containing a dictionary |
| `-c` | `--context` | How much context should be added at up and down stream of the sequence |
| `-n` | `--experiment_name` | Name of the data source of the RRIs e.g. experiment and organism |
| `-p` | `--param_file` | IntaRNA parameter file. Default: file in rrieval/IntaRNA_param |
| `-st` | `--use_structure` | Set 'off' if you want to disable structure, default 'on' |
| `-on` | `--out_name` | Name for the output directory instead of the current data |
| `-hf` | `--hand_feat` | If you want to start from hand curated feature files meaning starting already with positive and negative input files set 'on'. Default: 'off' |
| `-j` | `--n_jobs` | Number of jobs used for graph feature computation. Default: 1|


#### Output 
At the end of the run the location of the result table is given.
The final result table will have all columns of the input table and an additional prediction column, where you find the predicted class of each RRI (0 or 1).

Throughout the program several output files are generated and stored in the following structure:

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
Within CheRRI's **train** mode you can train your own model. 
This input data are the RRI interactions found by Direct Duplex Detection (DDD) methods. To extract this RRI interaction form DDD methods a tool named ChiRA is used to generate the 'ChiRA interaction summary' table. CheRRI expects as input the 'ChiRA interaction summary' file.

If you want to prepare a 'ChiRA interaction summary' table file please follow the [tutorial](https://training.galaxyproject.org/training-material//topics/transcriptomics/tutorials/rna-interactome/tutorial.html).

Starting from the RRI site information CheRRI will build a model based on features generated from the DDD method interactions site data. 

#### Retrieve RNA-RNA interactome files
Please select the ChiRA RRI output files as input for CheRRI **train** mode. Retive one file per replicate. You should specify the path to the folders containing all replicates. If you want to filter for interactions, where ChiRA already found a hybrid within its analysis you should enable the hybridization with the ChiRA workflow and set within CheRRI's  parameters filter_hybrid to 'on'. Some of the possible interaction maid be missed this way, but it will also filter out potential false positive interactions. 

#### Example call CheRRI train model mode
```
cherri train -i1 ./Cherri/test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7
```



#### Input Parameter
##### required
| ID | name | description |
|---|---|-----|
| `-i1` | `--RRI_path`| Path to folder storing the ChiRA interaction summary files|
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data, feature and model preparation |
| `-r` | `--list_of_replicates`| List the ChiRA interaction summary file for each replicate |
| `-l` | `--chrom_len_file`| Tabular file containing data in two columns format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-g` | `--genome`| Path to 2bit genome file, or used the build in download if you want the human or mouse genome |
##### optional
| ID | name | description |
|---|---|-----|
| `-c` | `--context`| How much context should be added at up and down stream of the sequence |
| `-n` | `--experiment_name`| Name of the data source of RRIs. Will be used for the file names|
| `-p` | `--param_file`| IntaRNA parameter file. Default: file in rrieval/IntaRNA_param |
| `-st` | `--use_structure`| Set 'off' if you want to disable graph-kernel features default: 'on' (when set to 'on' the feature optimization will be performed directly and the data will be stored in feature_files and no model/feature folder will be set up)|
| `-i2` | `--RBP_path`| Path to the genomic RBP cross link positions location (in BED format) |
| `-t` | `--run_time`| Time used for the optimization in seconds default: 12h|
| `-me` | `--memoryPerThread`| Memory in MB which each thread can use (total ram/threads)|
| `-j` | `--n_jobs`| Number of jobs used for graph feature computation and model selection. Default: 1|
| `-mi` | `--mixed`| Use mixed model to combine different dataset into a combined model. Default. 'off' (advanced option)| 
| `-fh` |`--filter_hybrid`| Filter the data for hybrids already detected by ChiRA (set 'on' to filter default:'off') |

For the mixed mode parameter you need to add all pre-computed models, who's data you would like to combine in one single folder (set -mi 'on'). Rename the first level folder to the experiment name and specify the path to this folder as '-i1' parameter and the experiment names as replicates (-r). The CheRRI pipeline will than first concatenate all positive and negative data set file of the given experiments and stars the CheRRI call form the feature generation step. 


#### Output 
At the end of the run the location of the trained model is given.

Throughout the program several output files are generated in the following structure:

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
Here we search for trusted RRIs, so RRIs which can be found in all replicates. In a first filter step only uniquely mapped RRIs are taken. Than RRI sequence partners in all replicates are found, using a overlap threshold. Output are the ChiRA input tables, now containing only the trusted RRIs. Out of all RRI pairs of the replicates only the one with the highest overlap to all others is added to the trusted_RRI data set. 

#### Input Parameter
| ID | name | description |
|---|---|-----|
| `-i` | `--input_path` | Path to folder storing input data (meaning storing all replicates) |
|`-r`| `--list_of_replicats` | List of filenames for all replicates |
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



