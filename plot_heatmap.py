#!/usr/bin/env python
import pandas as pd
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import rrieval.lib as rl
import re
import subprocess

import numpy as np

import seaborn as sns

def get_neg_data_name(full_name, data, context):
    prefix = data + '_context_method_'
    suffix = '_with_' + context + '_context_'
    name_temp = full_name.replace(suffix,'').replace(prefix,'')
    name_list = name_temp.split('_')
    name = name_list[0] + name_list[3]
    #name = re.sub(suffix, '', full_name)
    return name

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True
                                           , help= "path to input files")
    parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", required=True,
                        help= "path to output dir")
    parser.add_argument("-c", "--context", action="store", dest="context", required=True,
                        help= "context")
    parser.add_argument("-d", "--data_name", action="store", dest="data_name", required=True,
                        help= "name of dataset")



    args = parser.parse_args()

    input_file = args.input_file
    out_dir = args.out_dir
    context = args.context
    data_name = args.data_name

    df_input = pd.read_table(input_file, sep=',')

    model_list = list(df_input.columns.values)[4:]
    print(model_list)


    df_input['neg_data'] = df_input['suffeling_id'].apply(lambda x: get_neg_data_name(x, data_name, context))

    df_input['feature_id'] = 'nan'

    for idx, id in enumerate(df_input['featurs'].unique()):
        print('\nfeatures: %s\nid:%i\n'%(id, idx))
        df_input.loc[df_input.featurs == id, 'feature_id'] = idx

    #print(df_input['feature_id'])
    # index features
    # values AUC
    # colum datasets

    #print(df_input['neg_data'])
    COLUMN_NAMES = ['featurs', 'neg_data', 'AUC']

    count=0
    n=8

    fig, axes = plt.subplots(n, 1)
    fig.suptitle('Models AUC')

    for cname in model_list:
        temp_list = cname.split('_')
        #print(temp_list[0])
        if temp_list[1] == 'auc':
            count+=1
            model = temp_list[0]

            data = {'featurs' : df_input['feature_id'],
            'neg_data': df_input['neg_data'],
            'AUC': df_input[cname]}
            df_plot = pd.DataFrame(data, columns=COLUMN_NAMES)
            #print(df_plot)
            #file_name= out_dir + 'heatmap_' + model + '.pdf'

            #### plot
            #sns.set()
            #plt.subplot(n,1,count)
            #plt.figure(figsize=[5,5])
            df_plot = df_plot.pivot("featurs", "neg_data", "AUC")
            if count == n:
                ax = sns.heatmap(df_plot, ax=axes[(count-1)], vmin=0.4, vmax=1)
            else:
                ax = sns.heatmap(df_plot, xticklabels=False, ax=axes[(count-1)], vmin=0.4, vmax=1)

            #axes[(count-1)].set_title(model, fontsize='small')
            axes[(count-1)].tick_params(axis='both', labelsize=8)
            #axes[(count-1)].set_ticklabels([])
            #axes[(count-1)].set_yticklabels([])
            #axes[(count-1)].set_yticks([])
            axes[(count-1)].set_ylabel(model, fontsize=4, rotation=75)
            #plt.xticks(fontsize=6, rotation=90)




    file_name = out_dir +  data_name + '_heatmap.pdf'
    plt.savefig(file_name)


    df_input.to_csv(out_dir + data_name + '_modle_AUC_updated.csv', index=False)



if __name__ == '__main__':
    main()
