# Test Datasets

### Training data:
To thest if the Cherri **train** mode is running for you you can use one of the given datasets in the test data. 
The ./training/Paris/miRNA_human[1-3].tabular files are subsets of the original data and well suited to test the opation of the Cherri tool.

```
cherri train -i1 test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7
```

Check localy on denbi cloud:
cherri train -i1 ./Cherri/test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7

cherri train -i1 /vol/scratch/data_storage/data/RRIs/SPLASH_without_hybrids/ -r ChiRA_interaction_summary_hES_1_snRNAs.tabular ChiRA_interaction_summary_hES_2_snRNAs.tabular -g human -l human -o ./ -n splash_snRNA_human_test -c 50 -st on -t 600 -me 8000 -j 7


### Evaluation data:
The ./evaluate/test_evaluate_rri.cvs dataset is a artifical data set created to test Cherris **eval** function.

```
cherri eval -i1 test_data/evaluate/test_evaluate_rris.cvs -g human -l human -o ./ -n test_eval -c 150 -st on -m Cherri_models_data/Model_with_graph_features/PARIS_human/model/full_PARIS_human_context_150.model -mp Cherri_models_data/Model_with_graph_features/PARIS_human/feature_files/training_data_PARIS_human_context_150.npz -i2 Cherri_models_data/Model_with_graph_features/PARIS_human/occupied_regions/occupied_regions.obj
```



### Input dataset

#### full ChiRA extrat table header
```
'#reads','chrom_1st','start_1st','end_1st', 'strand_1st',
                'chrom_2end','start_2end','end_2end', 'strand_2end',
                'interaction_site_1st', 'interaction_site_2end',
                'IntaRNA_prediction', 'energy',
                'seq_1st_interaction_site', 'seq_2end_interaction_site',
                'start_interaction',
                'chrom_seq_1st_site', 'start_seq_1st_site',
                'stop_seq_1st_site','strand_seq_1st_site',
                'chrom_seq_2end_site', 'start_seq_2end_site',
                'stop_seq_2end_site','strand_seq_2end_site',
                'TPM_seq_1st_site', 'TPM_seq_2end_site', 'TPM_summary',
                'score_seq_1st_site', 'score_seq_2end_site','score_product',
                'biotype_region_1st', 'biotype_region_2end', 'ID_1st','ID_2end'
```

#### Header for user defined data 
```
'chrom_1st','start_1st','end_1st','strand_1st','chrom_2end','start_2end',
'end_2end','strand_2end'
```
The following coloums are needed within the cherri prediction but can be set withein the data construction. 
```
'IntaRNA_prediction','score_seq_1st_site','score_seq_2end_site'
```
#### If you with to reconstruct the CHirA input the following coloums also need to be present

Additional if a IntaRNA interaction for the interaciton sites is calculated pleas add

```
'interaction_site_1st', 'interaction_site_2end', 'IntaRNA_prediction', 'energy',
'seq_1st_interaction_site', 'seq_2end_interaction_site', 'start_interaction', 
'chrom_seq_1st_site', 'start_seq_1st_site', 'stop_seq_1st_site','strand_seq_1st_site',
'chrom_seq_2end_site', 'start_seq_2end_site', 'stop_seq_2end_site','strand_seq_2end_site'
```

Futher information about the interaction can be proviede in the following coloums
```
'TPM_seq_1st_site', 'TPM_seq_2end_site', 'TPM_summary', 'score_seq_1st_site', 
'score_seq_2end_site','score_product', 'biotype_region_1st', 'biotype_region_2end', 'ID_1st','ID_2end'
```


The file should be tabular separated!
