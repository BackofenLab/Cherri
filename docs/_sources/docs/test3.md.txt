# Documentation

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
