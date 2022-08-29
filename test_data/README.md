# Test Datasets

## Training data:
To thest if the Cherri **train** mode is running for you you cna use one of the given datasets in the test data. 
The ./training/Paris/miRNA_human[1-3].tabular are a subset of the original data and well suited to test the opation of the Cherri tool.

```
cherri train -i1 test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o ./ -n paris_human_test -c 50 -st on -t 600 -me 8000 -j 7
```


### Test data:
The ./evaluate/test_evaluate_rri.cvs dataset is a artifical data set created to test Cherris **eval** function.

```
cherri eval -i1 test_data/evaluate/test_evaluate_rris.cvs -g human -l human -o ./ -n test_eval -c 150 -st on -m Cherri_models_data/Model_with_graph_features/PARIS_human/model/full_PARIS_human_context_150.model -mp Cherri_models_data/Model_with_graph_features/PARIS_human/feature_files/training_data_PARIS_human_context_150.npz -i2 Cherri_models_data/Model_with_graph_features/PARIS_human/occupied_regions/occupied_regions.obj
```
