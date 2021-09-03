import pandas as pd
import math
#import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict
from interlap import InterLap
import sys
import argparse
import numpy as np
import pickle
import seaborn as sns

def append_to_list(in_list, append_list):
    for i in in_list:
        append_list.append(i)
    return append_list


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_path", action="store", dest="input_path",
                        required=True,
                        help= "path to folder storing all input data")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of positve trusted RRIs")



    args = parser.parse_args()
    input_path = args.input_path
    experiment_name = args.experiment_name

    overlap_thresholds = ['0.3', '0.4', '0.5', '0.6', '0.7']
    overlap_avg_val_list = []
    overlap_th_list = []
    len_avg_val_list = []


    for overlap_th in overlap_thresholds:
        #print(overlap_th)
        file_over = input_path + experiment_name + 'avg_overlap_' + overlap_th + '.obj'
        file_len = input_path + experiment_name + 'len_overlap_' + overlap_th + '.obj'
        #print(file)
        overlap_handle = open(file_over,'rb')
        overlap_avg_val = pickle.load(overlap_handle)
        # print(overlap_avg_val)
        overlap_handle.close()
        overlap_avg_val_list = append_to_list(overlap_avg_val, overlap_avg_val_list)
        overlap_th_temp = [overlap_th]*len(overlap_avg_val)
        overlap_th_list = append_to_list(overlap_th_temp, overlap_th_list)

        len_handle = open(file_len,'rb')
        len_avg_val = pickle.load(len_handle)
        # print(overlap_avg_val)
        len_handle.close()
        len_avg_val_list = append_to_list(len_avg_val, len_avg_val_list)

    #print(overlap_avg_val_list)
    #print(overlap_th_list)
    d = {'overlap_avg_val': overlap_avg_val_list, 'overlap_th': overlap_th_list, 'len_avg_val': len_avg_val_list}
    df_overlaps = pd.DataFrame(data=d)


    print(df_overlaps.info())


    # Plotting
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots()
    sns.boxplot(x="overlap_th", y="overlap_avg_val", data=df_overlaps, ax=ax)
    fig.savefig(input_path + '/avg_overlaps.pdf', format='pdf', dpi=300, bbox_inches='tight')

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots()
    sns.boxplot(x="overlap_th", y="len_avg_val", data=df_overlaps, ax=ax)
    fig.savefig(input_path + '/avg_len_overlaps.pdf', format='pdf', dpi=300, bbox_inches='tight')


    #### Plotting ######

    #histogrom enegy
    #fig1 = plt.figure()
    #bins = np.arange(min(enegy_list), max(enegy_list), 5)

    #plt.hist(enegy_list, bins=bins)
    #fig1.savefig(plot_path + "histogram_enegy.pdf", bbox_inches='tight')

    #seq1_len_list = [len[0] for len in interaction_length_also_nan]
    #seq2_len_list = [len[1] for len in interaction_length_also_nan]

    #d = {'rri_seq1': seq1_len_list, 'rri_seq2': seq2_len_list}
    #df_rri_len = pd.DataFrame(data=d)

    #myFig = plt.figure()
    #boxplot = df_rri_len.boxplot(column=['rri_seq1', 'rri_seq2'])

    #myFig.savefig(plot_path + "boxplot_rri_len_seq.pdf", bbox_inches='tight')

    # input_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/RNA_RNA_binding_evaluation/data/training/Paris/'
    # list_of_replicats = ['test_rep1.tabular', 'test_rep2.tabular', 'test_rep3.tabular']

if __name__ == '__main__':
    main()
