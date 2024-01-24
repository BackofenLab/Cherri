# CheRRI (Computational Help Evaluating RNA-RNA interactions)

CheRRI detects functional RNA-RNA interaction (RRI) sites, by evaluating if an interaction site most likely occurs in nature.
It helps to filter interaction sites generated either experimentally or by an RRI prediction algorithm by removing false positive interactions.

It is an open source project [hosted on GitHub](https://github.com/BackofenLab/Cherri).


## Table of content
1. [Installation](#Installation)
    - [Install Conda](#Install-Conda)
    - [Create environment manually](#Create-environment-manually)
    - [Manual installation](#Manual-installation)
    - [Install using pip](#Install-using-pip)
    - [Install CheRRI Conda package](#Install-CheRRI-Conda-package)
2. [Usage](#Usage)
    - [Evaluation of RRIs](#Evaluation-of-RRIs)
    * [RRI input format in evaluation mode](#RRI-input-format-in-evaluation-mode)
    * [Example call for CheRRI's evaluation mode](#Example-call-for-CheRRIs-evaluation-mode)
    * [Input parameters in evaluation mode](#Input-parameters-in-evaluation-mode)
    * [Output in evaluation mode](#Output-in-evaluation-mode)
    * [Validate your model using the **eval** mode](#Validate-your-model-using-the-eval-mode)
    - [Build a new CheRRI model in training mode](#Build-a-new-CheRRI-model-in-training-mode)
    * [Retrieve RNA-RNA interactome files using ChiRA](#Retrieve-RNA-RNA-interactome-files-using-ChiRA)
    * [Build RRI interactome file as input for CheRRI](#Build-RRI-interactome-file-as-input-for-CheRRI)
    * [Example call for CheRRI's training mode](#Example-call-for-CheRRIs-training-mode)
    * [Input parameters in training mode](#Input-parameters-in-training-mode)
    * [Run train in mixed model mode](#Run-train-in-mixed-model-mode)
3. [CheRRI's core method scripts](#cherris-core-method-scripts))
    - [RRI detection with find_trusted_RRI.py](#rri-detection-with-find_trusted_rripy))
    - [Compute occupied regions with find_occupied_regions.py](#compute-occupied-regions-with-find_occupied_regionspy))
    - [Interaction predictions with generate_pos_neg_with_context.py](#interaction-predictions-with-generate_pos_neg_with_contextpy))
    * [IntaRNA parameters used within CheRRI](#intarna-parameters-used-within-cherri))
    - [Feature extraction with get_features.py](#feature-extraction-with-get_featurespy))
    - [Feature selection and optimization](#feature-selection-and-optimization)
    

## CheRRI's workflow

CheRRI can be run in two modes, the model generation or **train** mode, or the RRI evaluation or **eval** mode. 
For the evaluation of a given set of RRI sites, a model must be specified in CheRRI's **eval** mode. Here pre-trained models can be applied or the user trains a model using CheRRI's **train** mode.
To train a novel model, an RNA-RNA interactome dataset specifying all RRI sites should be provided. CheRRI makes use of replicate data by checking if an RRI site can be found in all replicates within an overlap threshold. This is how CheRRI builds the set of trusted RRIs. In **eval** mode, the interaction positions are reformatted in the same way as the trusted RRIs. In both modes, CheRRI uses the same core method to generate a feature set which will be used to select a model in **train** mode, and in **eval** mode for the evaluation of the biological relevance of the submitted RRI sites.


[<img src="./plots/CheRRI-workflow.png" width="70%" />](./plots/CheRRI-workflow.png)



## Installation

CheRRI was developed in Linux and tested on Ubuntu (18.04 LTS). Conda is required to install CheRRI.


### Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the newest Python 3 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be accessed on the command line via (note that your version can be different):

```
$ conda --version
conda 4.10.3
```



### Create environment manually

To manually install CheRRI, first create a Conda environment:

```
conda create -n cherri python=3.8 -c conda-forge -c bioconda
conda activate cherri
```
Inside the environment, you need to install the following dependencies:


```
conda install -c conda-forge scikit-learn
conda install -c conda-forge networkx
conda install -c bioconda bedtools
conda install -c conda-forge biopython
conda install -c conda-forge interlap
conda install -c bioconda intarna
conda install -c conda-forge numpy
conda install -c conda-forge pandas
conda install -c conda-forge eden-kernel
conda install -c conda-forge biofilm
conda install -c conda-forge python-wget
```

Or create the environment with all dependencies at once:

```
conda create -n cherri -c conda-forge -c bioconda scikit-learn networkx numpy bedtools biopython interlap pandas intarna eden-kernel biofilm python-wget
conda activate cherri
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
conda activate cherri
```

### Manual installation

To install the tool itself, simply clone the repository and execute the installation script inside the cloned folder:

```
git clone https://github.com/BackofenLab/Cherri.git 
```

Make sure you are inside the CheRRI's conda environment and you are inside the tools folder. 

```
conda acivate cherri
cd Cherri
```

Than you can install CheRRI.

```
python -m pip install . --ignore-installed --no-deps -vv 
```


Now you can run CheRRI from any given folder:

```
cherri -h
```

### Install using pip
If you don't want to download the CheRRI's git folder you can also use the pipy package. 
```
pip install cherri
```

### Install CheRRI Conda package

work in progress


## Usage

You can use CheRRI in two modes. The **eval** mode predicts whether an RRI site of interest is biologically relevant. Pre-computed human and mouse models exist and can be downloaded from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.6533931). If you would like to build a model based on novel RRI interactome data, you can use the **train** mode.

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

We recommend that the user also specifies the occupied regions used to train the model provided via `--model_file`. For example, for the PARIS_human model without graph features, the occupied regions data object is located at 
`Model_without_graph_features/PARIS_human/occupied_regions/occupied_regions.obj` (inside the mentioned zip file from [Zenodo](https://doi.org/10.5281/zenodo.6533932)).




#### Example call for CheRRI's evaluation mode

For the test call please download the [Cherri_models_data](https://doi.org/10.5281/zenodo.6533932) zip folder. The PARIS_human model is needed to execute the call. Be sure to provide the correct location for the model and its feature set (`-m`, `-mp`). For example, assuming the data (zip folder extracted to folder `Cherri_models_data`) is stored inside the CheRRI folder:

```
cherri eval -i1 test_data/evaluate/test_evaluate_rris.csv -g human -l human -o ./ -n test_eval -c 150 -st on -m Cherri_models_data/Model_with_graph_features/PARIS_human/model/full_PARIS_human_context_150.model -mp Cherri_models_data/Model_with_graph_features/PARIS_human/feature_files/training_data_PARIS_human_context_150.npz -i2 Cherri_models_data/Model_with_graph_features/PARIS_human/occupied_regions/occupied_regions.obj
```



#### Input parameters in evaluation mode

Input parameters for CheRRI's **eval** mode (`cherri eval`):

##### Required:
| ID | name | description |
|---|---|-----|
| `-i1` |`--RRIs_table` | Table containing all RRIs that should be evaluated in the correct format|
| `-g` | `--genome_file`| Path to genome FASTA file, or use the built-in download function if you want the human or mouse genome |
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data and feature preparation as well as the evaluated instances |
| `-l` | `--chrom_len_file` | Tabular file containing data in two-column format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-m` | `--model_file` | Set path to the model which should be used for evaluation |
| `-mp` | `--model_params` | Set path to the feature file of the given model |
##### Optional:
| ID | name | description |
|---|---|-----|
| `-i2` | `--occupied_regions` | Path to occupied regions python object file containing a dictionary |
| `-c` | `--context` | How much context should be added at up- and downstream of the sequence. Default: 150 |
| `-n` | `--experiment_name` | Name of the data source of the RRIs, e.g. experiment and organism. Default: eval_rri |
| `-p` | `--param_file` | IntaRNA parameter file. Default: file in path_to_cherri_folder/Cherri/rrieval/IntaRNA_param |
| `-st` | `--use_structure` | Set 'off' if you want to disable structure. Default 'on' |
| `-on` | `--out_name` | Name for the output directory. Default: 'date_Cherri_evaluating_RRIs' |
| `-ef` | `--eval_features` | If you want to start from hand-curated feature files. Use this for evaluating test set performance (set 'on'). Default: 'off' |
| `-j` | `--n_jobs` | Number of jobs used for graph feature computation. Default: 1|


#### Output in evaluation mode

At the end of the run the location of the results table is given.
The final results table will have your the query and target ID's or your input sequences (`target_ID`,`query_ID`), the score of your instance (`instance_score`), the predicted class of each RRI (0 or 1) (`predicted_label`), if you are running the validation mode with `-ef on` the positive or negative label is given (`true_lable`), and finally all features of the instance are provided.

The Ids are a summary of `chromosme;strand;start;stop` oft the first (target) and the second (query) sequence.

Throughout the program, several output files are generated and stored in the following structure:

    ├── date_Cherri_evaluation_mode
    |   ├── evaluate_RRIs.csv
    |   ├── date_occ_out
    |       ├── occupied_regions.obj
    |       ├── rri_occupied_regions_overlapTH_0.3_scoreTH_1.csv
    |   ├── positive_instance
    |       ├── {name}_context_{context}pos.csv
    |       ├── {name}_context_{context}_block_ends_0_RRI_dataset.csv
    |   ├── feature_files
    |       ├── feature_filtered_{name}_context_{context}_pos.csv
    |       ├── training_data_{name}_context_{context}.npz
    |   ├── evaluation
    |       ├── evaluation_results_{name}.csv


#### Validate your model using the **eval** mode
You can also use CheRRIs **eval** mode to create a validation result table and than use the [compute_f1](./scripts/plots/compute_f1.py) to get the F1 score.

In the following is a example call to validate a theoretical model build from DataA with data form a different source e.g. DataB
```
cherri eval -i1 /path/to/Model_folder/Datab/feature_files/feature_filtered_<DataB>_context_<150>_pos_occ -g human -l human -o /path/to/Model_folder -n <eval_modelA_using_DataB> -c 150 -st on -m  /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model -mp  /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz -j 10 -on evaluation -ef on
```

In the following is a example call to cross validate a theoretical model build from DataA. The biofilm-cv will split the data into 5 parts used 4 to retrain a model and the left out to evaluate. 


```
python -m biofilm.biofilm-cv --infile /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 0 --model /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model --out /path/to/Model_folder
```

```
python -m biofilm.biofilm-cv --infile /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 1 --model /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model --out /path/to/Model_folder
```
```
python -m biofilm.biofilm-cv --infile /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 2 --model /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model --out /path/to/Model_folder
```
```
python -m biofilm.biofilm-cv --infile /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 3 --model /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model --out /path/to/Model_folder
```
```
python -m biofilm.biofilm-cv --infile /path/to/Model_folder/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 4 --model /path/to/Model_folder/DataA/model/optimized/full_<DataA>_context_<150>.model --out /path/to/Model_folder
```



Than use the result file to compute the F1 score using [compute_f1](./scripts/plots/compute_f1.py).

### Build a new CheRRI model in training mode

Within CheRRI's **train** mode you can train your own model. 
The input data are the RRIs found by Direct Duplex Detection (DDD) methods or other interactome HTS protocols. In theory it can be any RRI which is enough to build a solid model training dataset. If the interactions are from different organisms CheRRI needs to be called in a mixed_model 'on' training mode (will be explained later).

#### Retrieve RNA-RNA interactome files using ChiRA

To extract RRI interactions from DDD methods, a tool named [ChiRA](https://github.com/pavanvidem/chira) is used to generate the 'ChiRA interaction summary' table. CheRRI expects as input the 'ChiRA interaction summary' file.

If you want to prepare a 'ChiRA interaction summary' table file, please follow this [tutorial](https://training.galaxyproject.org/training-material//topics/transcriptomics/tutorials/rna-interactome/tutorial.html). You should prepare one ChiRA interaction summary file per replicate.

Starting from the RRI site information, CheRRI will build a model based on features generated from the DDD method interactions site data. 

ChiRA RRI output files are needed as input for CheRRI **train** mode. `--RRI_path` (`-i1`) demands the path to the the ChiRA interaction summary files, and `--list_of_replicates` (`-r`) demands the ChiRA interaction summary file names of the replicates used by CheRRI inside the `-i1` folder.

#### Build RRI interactome file as input for CheRRI
You can add your own interaction data in a tabular of csv format. Please follow the specific guides given below to build the table. 
You would need to provide the position information on an interaction using the following header line:
```
['chrom_1st','start_1st','end_1st','strand_1st','chrom_2end','start_2end','end_2end','strand_2end']
```

If you want to add the result of an interaction you should set all interaction position by using the following header names:
predicted interaction positions:
```
'chrom_seq_1st_site', 'start_seq_1st_site','stop_seq_1st_site','strand_seq_1st_site','chrom_seq_2end_site', 'start_seq_2end_site', 'stop_seq_2end_site','strand_seq_2end_site'
```
This information of the interaction need to be complete or not provided at all. If they are not provided the RRI position information from above are used as score for finding overlaps between replicates. 

If you have a score for the interactions you can also provide it in the following columns:
```
'score_seq_1st_site', 'score_seq_2end_site'
```

You can check the [example file](./test_data/training/user_defined.csv) to get a impression how it could look.

#### Example call for CheRRI's training mode

This is an example call to evoke CheRRI's model training mode inside the CheRRI folder:

```
cherri train -i1 test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7
```


#### Input parameters in training mode 


Input parameters for CheRRI's **train** mode (`cherri train`):

##### Required:
| ID | name | description |
|---|---|-----|
| `-i1` | `--RRI_path`| Path to folder storing the ChiRA interaction summary files|
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data, feature and model preparation |
| `-r` | `--list_of_replicates`| List the ChiRA interaction summary file for each replicate |
| `-l` | `--chrom_len_file`| Tabular file containing data in two-column format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-g` | `--genome`| Path to genome FASTA file, or use the built-in download function if you want the human or mouse genome |
##### Optional:
| ID | name | description |
|---|---|-----|
| `-c` | `--context`| How much context should be added at up- and downstream of the sequence. Default: 150 |
| `-n` | `--experiment_name`| Name of the data source of RRIs. Will be used for the file names. Default: 'model_rri'|
| `-p` | `--param_file`| IntaRNA parameter file. Default: file in path_to_cherri_folder/Cherri/rrieval/IntaRNA_param |
| `-st` | `--use_structure`| Set 'off' if you want to disable graph-kernel features. Default: 'on' (when set to 'on' the feature optimization will be performed directly and the data will be stored in feature_files and no model/feature folder will be created)|
| `-i2` | `--RBP_path`| Path to the genomic RBP crosslink or binding site locations (in BED format) |
| `-t` | `--run_time`| Time used for the optimization in seconds, default: 43200 (12h) |
| `-me` | `--memoryPerThread`| Memory in MB each thread can use (total ram/threads). Default: 4300|
| `-j` | `--n_jobs`| Number of jobs used for graph feature computation and model selection. Default: 1|
| `-mi` | `--mixed`| Use mixed model to combine different datasets into a combined model. Default: 'off' | 
| `-fh` |`--filter_hybrid`| Filter the data for hybrids already detected by ChiRA (set 'on' to filter). Default: 'off' |
| `-on` |`--out_name`| Name for the output directory. Default 'date_Cherri_evaluating_RRIs' |
| `-tp` |`--temp_dir`| Set a temporary directory for autosklearn. Either proved a path or 'out' to set it to the output directory. Default: 'off' |
| `-so` |`--no_sub_opt`| # of interactions IntraRNA will give is possible. Default: 5|
| `-es` |`--exp_score_th`|score threshold for the additional occupied regions [BED] Default: 10|
| `-ca` |`--context_additional`| context to extend left and right for the BED file instances. Default: 5|
| `-cv` |`--do_cv`| 5-fold cross validated of the pipeline will be performed using the training data. Set 'off' to skip. Default: 'on'|



#### Output in training mode

At the end of the run the location of the trained model is given.

Throughout the program, several output files are generated inside the output folder (default: `date_Cherri_model_build`), with the following structure:

    ├── date_CheRRI_model_build
    |   ├── date_occ_out
    |       ├── occupied_regions.obj
    |       ├── rri_occupied_regions_overlap_0.3
    |       ├── rri_occupied_regions_overlapTH_0.3_scoreTH_1.csv
    |   ├── pos_neg_data
    |       ├── {name}_context_{context}_pos_occ__block_ends_40_RRI_dataset.csv
    |       ├── {name}_context_{context}_pos_occ_neg.csv
    |       ├── {name}_context_{context}_pos_occ_pos.csv
    |   ├── feature_files
    |       ├── feature_filtered_{name}_context_{context}_pos.csv
    |       ├── feature_filtered_{name}_context_{context}_neg.csv
    |       ├── training_data_{name}_context_{context}.npz
    |   ├── model
    |       ├── {name}_context_{context}_fold{0-4}.model
    |       ├── {name}_context_{context}_fold{0-4}.csv
    |       ├── full_{name}_context_{context}.model


#### Run train in mixed model mode

CheRRI is able to build on model based on different datasets. Is the mixed parameter is set to 'on' CheRRI will connect training data for different datasets. However before running the mixed mode one would create the training data for the individual datasets.
Next we have a theoretical example of DataA, DataB and DataC, which should be trained together
Therefore best create an dedicated output folder e.g. CheRRI_build_model

Than build you first model:
```
cherri train -i1 /path/to/CheRRI_build_model/ -r dataA_1.tabular dataA_2.tabular -g human -l human -o ./ -n Data_A -c 150 -st on -t 600 -me 8000 -j 7
```
If you don't need the model of the individual datasets you can either set the -t very low or even interrupt the call once the model is build.
Next rename the output folder created by CheRRI to the name you gave to the Data/Model (-n Data_A)
```
mv /path/to/CheRRI_build_model/<date>_CheRRI_build_model /path/to/CheRRI_build_model/Data_A
```

Than run the next dataset:
```
cherri train -i1 /path/to/CheRRO_build_model/ -r dataB_1.tabular dataB_2.tabular -g human -l human -o ./ -n Data_B -c 150 -st on -t 600 -me 8000 -j 7
```
Next rename the output folder created by CheRRI to the name you gave to the Data/Model (-n Data_A)
```
mv /path/to/CheRRI_build_model/<date>_CheRRI_build_model /path/to/CheRRI_build_model/Data_B
```

Than run the last dataset:
```
cherri train -i1 /path/to/CheRRI_build_model/ -r dataC_1.tabular dataC_2.tabular -g mouse -l mouse -o ./ -n Data_C -c 150 -st on -t 600 -me 8000 -j 7
```
And rename the output folder created by CheRRI to the name you gave to the Data/Model (-n Data_A)
```
mv /path/to/CheRRI_build_model/<date>_CheRRI_build_model /path/to/CheRRI_build_model/Data_C
```

Finally you can run CheRRI **train** in the mixed model mode like this:
```
cherri train -i1 /path/to/CheRRI_build_model/  -r Data_A Data_B Data_C  -g /not/needed/ -l /not/needed/ -o /path/to/CheRRI_build_model/ -n Full 
```
This time your replicates are the names of the training datasets you want to connect (Data_A Data_B Data_C). 



## CheRRI's core method scripts

CheRRI is built as a modular tool calling individual scripts accomplishing the various tasks of CheRRI's core methods. If you would like to perform only one step of the CheRRI pipeline, you can do this by calling the individual scripts. A short description of this scripts is given in the following.


### RRI detection with find_trusted_RRI.py
Here we search for trusted RRIs, so RRIs which can be found in all replicates. In a first filter step only uniquely mapped RRIs are taken. Than RRI sequence partners in all replicates are found, using a overlap threshold. Output are the ChiRA input tables, now containing only the trusted RRIs. Out of all RRI pairs of the replicates only the one with the highest overlap to all others is added to the trusted_RRI data set. 

#### Input parameters for find_trusted_RRI.py
| ID | name | description |
|---|---|-----|
| `-i` | `--input_path` | Path to folder storing input data (containing all replicates) |
|`-r`| `--list_of_replicats` | List of file names for all replicates |
| `-o` | `--overlap_th` | Overlap threshold to find trusted RRIs |
| `-d` | `--output_path` | Path where output folder should be stored |
|`-n` | `--experiment_name` | Name of the data source of positive trusted RRIs |
| `-s` | `--score_th` | Threshold for EM score from ChiRA |
| `-fh` | `--filter_hybrid` | Filter the data for hyprids alrady detected by ChiRA |

#### Output of find_trusted_RRI.py
The filtered set of trusted RRI sites in tabular format. 


### Compute occupied regions with find_occupied_regions.py
Given the RRI information tables from ChiRA and RNA-protein binding positions, an InterLab object is build. The occupied information can be used to mask parts of the genome and therefore enable to select negative interaction regions.


#### Input parameters for find_occupied_regions.py
| ID | name | description |
|---|---|-----|
|`-i1` | `--RRI_path` | Path to folder storing all RRI data (table) |
| `-i2` | `--rbp_path` | Path to RBP site data file (BED format) |
| `-r` | `--list_of_replicates` | List of file names for all replicates |
| `-o` | `--out_path` | Path where output folder should be stored |
| `-t` | `--overlap_th` | Overlap threshold |
| `-s` | `--score_th` | Score threshold | 
| `-e` | `--external_object` | External RRI overlapping object (InterLap dict)|
| `-fh` | `--filter_hybrind` | Filter the data for hybrids already detected by ChiRA |
| `-mo` | `--mode` | Function call within which CheRRI mode (train/eval)|

#### Output of find_occupied_regions.py

A python pickle file object storing occupied regions in an InterLap dictionary. 


### Interaction predictions with generate_pos_neg_with_context.py

Given a set of trusted RRI sites and occupied regions, a given context is appended. Then positive interactions are computed by calling IntaRNA, specifying the trusted RRI sites as seed regions. The negative interactions are computed by calling IntaRNA on regions outside the RRI sites / occupied regions.


#### Input parameters for generate_pos_neg_with_context.py

| ID | name | description |
|---|---|-----|
| `-i1` | `--input_rris` |Path to file storing all trusted RRIs|
|  `-i2` | `--input_occupied` | Path to occupied regions file |
| `-d` | `--output_path` | Path where output folder should be stored |
| `-n` | `--experiment_name` | Name of the data source of positive trusted RRIs|
| `-g` | `--genome_file` | Path to genome FASTA file |
| `-c` | `--context` | How much context should be added up- and downstream |
|   | `--pos_occ`  | Occupied regions are set (default) |
|   | `--no_pos_occ`  | Set if no occupied regions should be used |
| `-b` | `--block_ends` | # nucleotides blocked at the ends of each extended RRI site |
| `-s` | `--no_sub_opt` | # of interactions IntraRNA will give if possible |
| `-l` | `--chrom_len_file` | Tabular file containing chrom name \t chrom length for each chromosome |
| `-p` | `--param_file` | IntaRNA parameter file |
| `-m` | `--mode` | Which CheRRI mode is running (train/eval) |



#### IntaRNA parameters used within CheRRI

To generate the current features IntaRNA parameters by default are set to:

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


IntaRNA parameters can be changed by specifying a custom IntaRNA parameter file.



#### Output of generate_pos_neg_with_context.py
Positive and negative datasets stored in tabular format. 


### Feature extraction with get_features.py

Here additional sequence features are computed and the output is filtered for the final or a given feature set. The features are stored in tabular format.
Note that the graph-kernel features are computed directly in CheRRI's main functions and do not have a separate file. 

#### Input parameters for get_features.py

| ID | name | description |
|---|---|-----|
| `-i` | `--input` | Path to input file |
| `-f` | `--feature_set_list` | Set of features the script will output |
| `-o` | `--output_file` | Output file path including the file name |


#### Output of get_features.py

Tabular file having all features given via `--feature_set_list`.


### Feature selection and optimization

If you would only like to run the feature selection or optimization, please check out [biofilm](https://github.com/smautner/biofilm). 



