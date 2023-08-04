#!/bin/sh

### Input variabels

OUT="/vol/scratch/data_storage/test_cherri"
GENOMEFA="/paris_human_test1/genome/hg38.fa"
GENOMETABLE="/paris_human_test1/genome/hg38.chrom.sizes"

# set up test dir
echo $OUT
mkdir -p $OUT

### train

# train fist test set
cherri train -i1 ./test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular miRNA_human_3.tabular -g human -l human -o $OUT -n paris_human_test1 -c 50 -st on -t 600 -me 8000 -j 7 -on paris_human_test1 -tp out -fo 3

# train second test set
cherri train -i1 ./test_data/training/Paris/ -r miRNA_human_1.tabular miRNA_human_2.tabular -g $OUT$GENOMEFA -l $OUT$GENOMETABLE -o $OUT -n paris_human_test2 -c 50 -st on -t 500 -me 8000 -j 7 -on paris_human_test2 -tp out -cv off

# train third test set
cherri train -i1 ./test_data/training/Paris/ -r miRNA_human_3.tabular -g $OUT$GENOMEFA -l $OUT$GENOMETABLE -o $OUT -n paris_human_test3 -c 50 -st on -t 500 -me 8000 -j 7 -on paris_human_test3 -tp out -cv off


# train combinde model
cherri train -i1 $OUT  -r paris_human_test1 paris_human_test2 paris_human_test3 -g /not/needed/ -l /not/needed/ -o $OUT -n Full_model_test -c 50 -st on -mi on -t 500 -me 8000 -j 7 


### eval 

# model evaluation - test1 with test2
cherri eval -i1 $OUT/paris_human_test2/feature_files/feature_filtered_paris_human_test2_context_50_pos_occ -g $OUT$GENOMEFA -l $OUT$GENOMETABLE -o $OUT -n test_cross_model -c 50 -st on -m  $OUT/paris_human_test1/model/optimized/full_paris_human_test1_context_50.model -mp  $OUT/paris_human_test1/feature_files/training_data_paris_human_test1_context_50.npz -j 10 -on evaluation -ef on

# model cross validation for fold 3
python -m biofilm.biofilm-cv --infile $OUT/paris_human_test1/feature_files/training_data_paris_human_test1_context_50.npz --foldselect 3 --model $OUT/paris_human_test1/model/optimized/paris_human_test1_context_50.model --out $OUT/evaluation/evaluation/paris_human_test1_paris_human_test1_fold3.csv


# evalutation using test eval test data
cherri eval -i1 ./test_data/evaluate/test_evaluate_rris.csv  -g $OUT$GENOMEFA -l $OUT$GENOMETABLE -o $OUT -n test_eval -c 50 -st on -m $OUT/paris_human_test1/model/optimized/full_paris_human_test1_context_50.model -mp  $OUT/paris_human_test1/feature_files/training_data_paris_human_test1_context_50.npz 

