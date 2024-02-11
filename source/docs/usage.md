# Usage

You can use CheRRI in two modes. The **eval** mode predicts whether an RRI site of interest is biologically relevant. Pre-computed human and mouse models exist and can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.10555733). If you would like to build a model based on novel RRI interactome data, you can use the **train** mode.

## Evaluation of RRIs

Based on a tabular file containing chromosomal position data of the RRIs, CheRRI classifies if the interaction region is likely to be a biologically relevant one.

For the **eval** mode, a model and the filtered feature set have to be specified.
CheRRI has pre-trained models for human and mouse, which can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.10555733) (see `content.txt` inside the zip folder for details).
If there exists an RNA-RNA-interactome dataset for your preferred organism, we recommend to train your own organism-specific model using CheRRI's **train** mode. After training, the model can be than used in **eval** mode for the classification of your predicted RRI positions.


### RRI input format in evaluation mode

The RRI instances to be evaluated need to be given in CSV format (parameter `--RRIs_table`). The table needs the following header line:

```
chrom1,start1,stop1,strand1,chrom2,start2,stop2,strand2
```

Following the header line, each subsequent line represents an RRI, with chromosome ID (format: 1,2,3 ...), interaction start, interaction end, and strand ("+" or "-") of the two interacting partners. For example, you might want to evaluate the following three RRI sites:

```
19,18307518,18307539,-,14,90454500,90454521,+
X,109054541,109054590,+,9,89178539,89178562,-
10,123136102,123136122,+,5,1245880,1245902,+
```





### Example call for CheRRI's evaluation mode

For the test call please download the [Cherri_models_data](https://doi.org/10.5281/zenodo.10555733) zip folder. The human model is needed to execute the call. Be sure to provide the correct location for the model and its feature set (`-m`, `-mp`). For example, assuming the data (zip folder extracted to folder `Cherri_models_data`) is stored inside the CheRRI folder:

```
cherri eval -i1 test_data/evaluate/test_evaluate_rris.csv -g human -l human -o ./ -n test_eval -c 150 -st on -m Cherri_models_data/Model_with_graph_features/human/full_human_context_150.model -mp Cherri_models_data/Model_with_graph_features/human/training_data_human_context_150.npz 
```

If you would like to preform a test evaluation you can use the [test_file](https://github.com/BackofenLab/Cherri/blob/master/test_data/evaluate/test_evaluate_rris.csv) stored in CheRRIs github repository. 

### Input parameters in evaluation mode

Input parameters for CheRRI's **eval** mode (`cherri eval`):

#### Required:
| ID | name | description |
|---|---|-----|
| `-i1` |`--RRIs_table` | Table containing all RRIs that should be evaluated in the correct format|
| `-g` | `--genome_file`| Path to genome FASTA file, or use the built-in download function if you want the human or mouse genome |
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data and feature preparation as well as the evaluated instances |
| `-l` | `--chrom_len_file` | Tabular file containing data in two-column format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-m` | `--model_file` | Set path to the model which should be used for evaluation |
| `-mp` | `--model_params` | Set path to the feature file of the given model |
#### Optional:
| ID | name | description |
|---|---|-----|
| `-i2` | `--occupied_regions` | Path to occupied regions python object. This file should be used if there are regions which that should be blocked from interactions. One can create this file with the find_occupied_regions.py |
| `-c` | `--context` | How much context should be added at up- and downstream of the sequence. Default: 150 |
| `-n` | `--experiment_name` | Name of the data source of the RRIs, e.g. experiment and organism. Default: eval_rri |
| `-p` | `--param_file` | IntaRNA parameter file. Default: file in path_to_cherri_folder/Cherri/rrieval/IntaRNA_param |
| `-st` | `--use_structure` | Set 'off' if you want to disable structure. Default 'on' |
| `-on` | `--out_name` | Name for the output directory. Default: 'date_Cherri_evaluating_RRIs' |
| `-ef` | `--eval_features` | If you want to start from hand-curated feature files. Use this for evaluating test set performance (set 'on'). Default: 'off' |
| `-j` | `--n_jobs` | Number of jobs used for graph feature computation. Default: 1|


### Output in evaluation mode

At the end of the run the location of the results table is given.
The final results table will have a query and target ID of your input sequences (`target_ID`,`query_ID`), the score of your instance (`instance_score`), the predicted class of each RRI (0 or 1) (`predicted_label`), if you are running the validation mode with `-ef on` the positive or negative label (`true_lable`), and subsequent all features of the instance are provided.

The IDs are a summary of `chromosome;strand;start;stop` oft the first (target) and the second (query) sequence.

Throughout the program, several output files are generated and stored in the following structure:

    ├── date_Cherri_evaluation_mode
    |   ├── evaluate_RRIs.csv
    |   ├── positive_instance
    |       ├── {name}_context_{context}pos.csv
    |       ├── {name}_context_{context}_block_ends_0_RRI_dataset.csv
    |   ├── feature_files
    |       ├── feature_filtered_{name}_context_{context}_pos.csv
    |       ├── training_data_{name}_context_{context}.npz
    |   ├── evaluation
    |       ├── evaluation_results_{name}.csv


### Validate your model using the **eval** mode
You can also use CheRRIs **eval** mode to create a validation result table and than use the [compute_f1](https://github.com/BackofenLab/Cherri/blob/master/scripts/plots/compute_f1.py) to get the F1 score.

In the following a example call to validate a theoretical model build from DataA with data form a different source e.g. DataB is given.

You can set a MODELPATH leading to your models e.g. of DataA and DataB. Here we assume that DataA and DataB will in the same MODELPATH location (e.g. `MODELPATH=~/myproject/cherri_models`).
You can set the MODELPATH on a linux based system like this:
```
export MODELPATH=$PATH:~/myproject/cherri_models
```

```
cherri eval -i1 $MODELPATH/<DataB>/feature_files/feature_filtered_<DataB>_context_<150>_pos_occ -g human -l human -o $MODELPATH -n <eval_modelA_using_DataB> -c 150 -st on -m  $MODELPATH/DataA/model/full_<DataA>_context_<150>.model -mp  $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz -j 10 -on evaluation -ef on
```

In the following is a example call to cross validate a theoretical model build from DataA. The biofilm-cv will split the data into 5 parts, use 4 to retrain a model and the 5th, left out part to evaluate the trained model. 


```
python -m biofilm.biofilm-cv --infile $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 0 --model $MODELPATH/DataA/model/optimized/full_<DataA>_context_<150>.model --out $MODELPATH
```

```
python -m biofilm.biofilm-cv --infile $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 1 --model $MODELPATH/DataA/model/optimized/full_<DataA>_context_<150>.model --out $MODELPATH
```
```
python -m biofilm.biofilm-cv --infile $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 2 --model $MODELPATH/DataA/model/optimized/full_<DataA>_context_<150>.model --out $MODELPATH
```
```
python -m biofilm.biofilm-cv --infile $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 3 --model $MODELPATH/DataA/model/optimized/full_<DataA>_context_<150>.model --out $MODELPATH
```
```
python -m biofilm.biofilm-cv --infile $MODELPATH/DataA/feature_files/training_data_<DataA>_context_<150>.npz --foldselect 4 --model $MODELPATH/DataA/model/optimized/full_<DataA>_context_<150>.model --out $MODELPATH
```



Than use the result file to compute the F1 score using [compute_f1](https://github.com/BackofenLab/Cherri/blob/master/scripts/plots/compute_f1.py).

## Build a new CheRRI model in training mode

Within CheRRI's **train** mode you can train your own model. 
The input data are the RRIs found by Direct Duplex Detection (DDD) methods or other interactome HTS protocols. In theory it can be any RRI which is enough to build a solid model training dataset. If the interactions are from different organisms CheRRI needs to be called in a mixed_model 'on' training mode (will be explained later).

### Retrieve RNA-RNA interactome files using ChiRA

To extract RRI interactions from DDD methods, a tool named [ChiRA](https://github.com/pavanvidem/chira) is used to generate the 'ChiRA interaction summary' table. CheRRI expects as input the 'ChiRA interaction summary' file.

If you want to prepare a 'ChiRA interaction summary' table file, please follow this [tutorial](https://training.galaxyproject.org/training-material//topics/transcriptomics/tutorials/rna-interactome/tutorial.html). You should prepare one ChiRA interaction summary file per replicate.

Starting from the RRI site information, CheRRI will build a model based on features generated from the DDD method interactions site data. 

ChiRA RRI output files are needed as input for CheRRI **train** mode. `--RRI_path` (`-i1`) demands the path to the the ChiRA interaction summary files, and `--list_of_replicates` (`-r`) demands the ChiRA interaction summary file names of the replicates used by CheRRI inside the `-i1` folder.

### Build RRI interactome file as input for CheRRI
You can add your own interaction data in a tabular CSV format. Please follow the specific guides given below to build the table. 
You would need to provide the position information on an interaction using the following header line:
```
['chrom_1st','start_1st','end_1st','strand_1st','chrom_2end','start_2end','end_2end','strand_2end']
```

If you want to add the result of an interaction you should set all interaction positions by using the following header names:
predicted interaction positions:
```
'chrom_seq_1st_site', 'start_seq_1st_site','stop_seq_1st_site','strand_seq_1st_site','chrom_seq_2end_site', 'start_seq_2end_site', 'stop_seq_2end_site','strand_seq_2end_site'
```
The information of the interaction need to be complete or not provided at all. If they are not provided the RRI position information from above are used as score for finding overlaps between replicates. 

If you have a score for the interactions you can also provide it in the following columns:
```
'score_seq_1st_site', 'score_seq_2end_site'
```

You can check the [example file](https://github.com/BackofenLab/Cherri/blob/master/test_data/training/user_defined.csv) to get a impression how it could look.

### Example call for CheRRI's training mode

This is an example call to evoke CheRRI's model training mode. If you download the Cherri github folder you can use the files for a test call:

```
cherri train -i1 test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n human_test -c 50 -st on -t 600 -me 8000 -j 7
```


### Input parameters in training mode 


Input parameters for CheRRI's **train** mode (`cherri train`):

#### Required:
| ID | name | description |
|---|---|-----|
| `-i1` | `--RRI_path`| Path to folder storing the ChiRA interaction summary files|
| `-o` | `--out_path`| Path to output directory where the output folder will be stored. It will contain separate output folders for each step of the data, feature and model preparation |
| `-r` | `--list_of_replicates`| List the ChiRA interaction summary file for each replicate |
| `-l` | `--chrom_len_file`| Tabular file containing data in two-column format for each chromosome: 'chrom name' \t 'chrom length'. You can directly specify 'human' or 'mouse' |
| `-g` | `--genome`| Path to genome FASTA file, or use the built-in download function if you want the human or mouse genome |
#### Optional:
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
| `-fo` |`--folds`| number of folds for the cross validation. Default: 5|
| `-mt` |`--methods`| Methods used for model selection. Default: any|



### Output in training mode

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


### Run train in mixed model mode

CheRRI is able to build one model based on different datasets. Is the mixed parameter is set to 'on' CheRRI will connect training data for different datasets. However before running the mixed mode one would create the training data for the individual datasets.
Next we have a theoretical example of DataA, DataB and DataC, which should be trained together
Therefore best create an dedicated output folder e.g. CheRRI_build_model

Then build you first model:
```
cherri train -i1 /path/to/CheRRI_build_model/ -r dataA_1.tabular dataA_2.tabular -g human -l human -o ./ -n Data_A -c 150 -st on -t 600 -me 8000 -j 7
```
If you don't need the model of the individual datasets you can either set the -t very low or even interrupt the call once the model is build.
Next rename the output folder created by CheRRI to the name you gave to the Data/Model (-n Data_A)
```
mv /path/to/CheRRI_build_model/<date>_CheRRI_build_model /path/to/CheRRI_build_model/Data_A
```

Then run the next dataset:
```
cherri train -i1 /path/to/CheRRO_build_model/ -r dataB_1.tabular dataB_2.tabular -g human -l human -o ./ -n Data_B -c 150 -st on -t 600 -me 8000 -j 7
```
Next rename the output folder created by CheRRI to the name you gave to the Data/Model (-n Data_A)
```
mv /path/to/CheRRI_build_model/<date>_CheRRI_build_model /path/to/CheRRI_build_model/Data_B
```

Then run the last dataset:
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




