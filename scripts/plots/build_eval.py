#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
import rrieval.lib as rl
from biofilm.util.data import loadfolds
from sklearn.metrics import f1_score


def call_script(call,reprot_stdout=False):
    """
    Starts a subprosses to call a script and checks for errors.


        Parameters
        ----------
        call : cmd comand

        Returns
        -------
        out
            if reprot_stdout set True returns the stdout (default: False)

    """
    process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
    #process.wait()
    out, err = process.communicate()
    #print(err.decode('utf-8'))
    error = err.decode('utf-8')

    #assert not error, "script is complaining:\n%s\n%s" %(call, error)
    if reprot_stdout == True:
        # out = out.decode('utf-8')
        return out

def remove_chr(chr):
    chr.replace('chr','')
    return chr


def create_eval_data(file):
    df_full = pd.read_csv(file)
    #print(pos_file)
    #print(df_pos.info())

    #df_pos['label'] = 1
    #df_neg['label'] = 0
    #df_full = pd.concat([df_pos,df_neg])

    header_eval = ['chrom1','start1','stop1','strand1','chrom2','start2','stop2','strand2']

    df_full['chrom1'] = df_full['target_key'].apply(lambda x: x.split(';')[0])
    df_full['start1'] =  (df_full['target_con_s'] + df_full['start1']).astype(int)
    df_full['stop1'] =  (df_full['target_con_s'] +  df_full['end1']).astype(int)
    df_full['strand1'] = df_full['target_key'].apply(lambda x: x.split(';')[1])
    df_full['chrom2'] = df_full['side_target'] = df_full['query_key'].apply(lambda x: x.split(';')[0])
    df_full['start2'] =  (df_full['query_con_s'] +  df_full['start2']).astype(int)
    df_full['stop2'] =  (df_full['query_con_s'] +  df_full['end2']).astype(int)
    df_full['strand2'] = df_full['query_key'].apply(lambda x: x.split(';')[1])

    df_full['chrom1'] = (df_full['chrom1'].apply(lambda x: x.replace('chr',''))).astype(str)
    df_full['chrom2'] = (df_full['chrom2'].apply(lambda x: x.replace('chr',''))).astype(str)

    return df_full[header_eval]

def get_eval_df(pos_file,neg_file):
    df_pos = pd.read_csv(pos_file)
    df_neg = pd.read_csv(neg_file)

    df_pos['label'] = 1
    df_neg['label'] = 0
    df_full = pd.concat([df_pos,df_neg])
    return df_full


def get_pos_neg_file_names(data_path,dataset):
    pos_file =  f'{data_path}pos_neg_data/{dataset}_context_150_pos_occ_pos.csv'
    neg_file =  f'{data_path}pos_neg_data/{dataset}_context_150_pos_occ_neg.csv'
    return pos_file, neg_file

def get_genome(organism):
    if organism == 'human':
        genome = '/vol/scratch/data_storage/data/genomes/hg38_UCSC_20210318.2bit'
        chrom_len = '/vol/scratch/data_storage/data/genomes/hg38_Info.tab'
    elif organism == 'mouse':
        genome = '/vol/scratch/data_storage/data/genomes/mm10_UCSC_20210324.2bit'
        chrom_len = '/vol/scratch/data_storage/data/genomes/mm10.chrom.sizes'
    return genome, chrom_len

def call_eval_cherri(file,label,genome, chrom_len,dataset,out_dir,in_path):
    # generate eval file
    out_eval_file = f'{out_dir}/{dataset}_{label}_eval.csv'
    df_eval = create_eval_data(file)
    df_eval.to_csv(out_eval_file, index=False)
    # call Cherri eval?
    call_eval = f'cherri eval -i1 {out_eval_file} -g {genome} -l {chrom_len} -o {out_dir} -n {dataset} -c 150 --use_structure on -mp {in_path}/{dataset}/model/features/{dataset}_context_150.npz -m {in_path}/{dataset}/model/optimized/{dataset}_context_150.model -on {dataset}'
    print(call_eval)
    out = call_script(call_eval,True)
    eval_file = get_filepath(out)
    return eval_file

def get_filepath(out):
    m = re.search('Result file: (.+?)evaluation_results', out)
    if m:
        found = m.group(1)
    #x.split('AAA')[1].split('ZZZ')[0]
    file = found + 'evaluation_results'
    return file

def calculate_measures(df_model):
    #df_model.columns = df_model.columns.str.strip()
    #print(df_model.columns)
    #print(df_model['true_label'])
    f1 = f1_score(df_model['label'].tolist(), df_model['prediction'].tolist(), average='macro')
    #print(f1)
    #precision, recall, thresholds, auc_prc = compute_prc(df_model['true_label'].tolist(), df_model['instance_score'].tolist())
    return f1


def filter_features(X,featurefile):
    ft = np.load(featurefile)['arr_0']
    #print(X.info())
    header = X.columns
    header_list = header.tolist()
    features_filterd=list(compress(header_list, ft))
    #print(features_filterd)
    X_filterd = X[features_filterd]
    #print('dfInfo after:')
    #print(X_filterd.info())

    return X_filterd


def classify(df_eval,in_model_filepath):
    #X = pd.read_csv(in_data_filepath, sep=',')
    #model_handle = open(in_model_filepath,'rb')
    #model = pickle.load(model_handle)
    #model_handle.close()
    #params = loadfile(in_model_filepath)['params']
    #print(params)
    model = loadfile(in_model_filepath)['estimator']
    y_pred=model.predict(df_eval)
    #print('model predictions')
    xtra = pd.DataFrame({'prediction': y_pred})
    df_result = pd.concat([df_eval, xtra], axis=1)
    # df_eval['prediction'] = y_pred
    #print(y_pred)

    return df_result


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-o", "--output_filepath", required=True, help= "output file")
    args = parser.parse_args()

    out_dir = args.output_filepath



    in_path = '/vol/scratch/data_storage/Cherri_build_model'
    # model type:
    model_pathes = ['Full', 'Full_human','PARIS_human','PARIS_human_RBPs','PARIS_mouse','SPLASH_human', 'Full_human_RRIs']
    #model_pathes = ['PARIS_human','PARIS_mouse']

    #organism_list = ['human','mouse']
    #f1_score = []

    for dataset_model in model_pathes:
        if dataset_model == 'PARIS_human':
            print('already computed!!!!')
            break
        X_model = f'{in_path}/{dataset_model}/feature_files/training_data_{dataset_model}_context_150.npz'
        model_feature = f'{in_path}/{dataset_model}/model/features/{dataset_model}_context_150'
        model_file = f'{in_path}/{dataset_model}/model/optimized/{dataset_model}_context_150.model'
        print(f'model set:{dataset_model}')
        #print(f'dataset model:{X_model}')
        #print(f'selected features model:{model_feature}')
        #print(f'model:{model_file}')
        call_refit = f'$(which python) -m biofilm.biofilm-cv --infile {X_model} --folds 0 --featurefile {model_feature} --model {model_file} --out {out_dir}/{dataset_model}'
        print(f'\n refit model:')
        print(call_refit)
        #call_script(call_refit)

        for dataset in model_pathes:
            print(f'data set:{dataset}')
            X_data = f'{in_path}/{dataset}/feature_files/training_data_{dataset}_context_150.npz'
            #print(f'\n test dataset:{X_data}')

            # X_filterd = filter_features(X,model_feature)
            if dataset_model == dataset:
                # do crossval
                print(f'\ncross val:\n')
                for i in range(5):

                    call_model = f"$(which python) -m biofilm.biofilm-cv --infile {X_data} --foldselect {i} --featurefile {model_feature} --model {model_file} --out {out_dir}/{dataset_model}_{dataset}_fold{i}"
                    print(call_model)
                    #call_script(call_model)

            else:
                #do cross model
                # refit
                print(f'\ncross model:\n')
                call_predict = f'$(which python) -m biofilm.util.out --folds 0 --featurefile {model_feature} --model {out_dir}/{dataset_model}.model --out {out_dir}/{dataset_model}_{dataset} --infile {X_data}'
                print(call_predict)
                #call_script(call_predict)
            #data_path = f'{in_path}/{dataset}/'

            #genome, chrom_len = get_genome(organism_list[idx])
            #pos_file, neg_file = get_pos_neg_file_names(data_path,dataset)

            #pos_eval_file = call_eval_cherri(pos_file,'pos',genome, chrom_len,dataset,out_dir,in_path)
            #neg_eval_file = call_eval_cherri(neg_file,'neg',genome, chrom_len,dataset,out_dir,in_path)

            #df_eval = get_eval_df(pos_eval_file,neg_eval_file)


    #(X,Y,_,_),_,_ = bf.loadfolds(X_file, folds=0, featurefile=model_feature) # somhow get data from here npz loadfolds (biofilm)!


            #df_eval = classify(X_filterd, model_file)





        #f1 = calculate_measures(df_eval)
        #f1_score.append(f1_score)

        #copy eval

    print('Results')
    #print(model_pathes)
    #print(f1_score)



        #out = '\n######################\nResult file: /some/path/fiel_evaluation_results\n############\n'
        #neg_eval_file = get_filepath(out)
        #print(neg_eval_file)


if __name__ == '__main__':
    main()
