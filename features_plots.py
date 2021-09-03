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


def plot_historgramm(colum, file_name, plot_dir, df_neg,df_pos):
    fig = plt.figure()
    n, bins, edges = plt.hist([df_neg[colum], df_pos[colum]], color=['r','b'], bins=15, alpha=0.5)
    plt.xticks(bins)
    plt.xticks(rotation=75)
    #print(min(df_neg[colum]))
    #print(max([-3,-5]))
    max_val = max([max(df_neg[colum]),max(df_pos[colum])])
    min_val = min([min(df_neg[colum]),min(df_pos[colum])])
    if min_val<0 and max_val <0:
        max_val_temp = max_val
        max_val = min_val
        min_val = max_val_temp
    plt.xlim(min_val, max_val)
    labels= ["negative", "positive"]
    plt.legend(labels)
    fig.savefig(plot_dir + file_name + '_histogramm.pdf', format='pdf', dpi=300, bbox_inches='tight')


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--pos_file",
                        help= "file to features of positive data",
                        default="/vol/scratch/data/features_files/test/fearutes_pos.csv")
    parser.add_argument("-i2", "--neg_file",
                        help= "file to features of negative data",
                        default="/vol/scratch/data/features_files/test/fearutes_neg.csv")



    args = parser.parse_args()

    pos_file = args.pos_file
    neg_file = args.neg_file
    file_add_name = 'True_b40_'

    plot_dir = "/vol/scratch/data/features_files/full_c300/" + file_add_name


    df_pos = pd.read_table(pos_file, sep=',')
    df_neg = pd.read_table(neg_file, sep=',')

    #print(df_pos.info())
    #print(df_neg.info())

    # print E, E_hybrid, ED1, ED2, no_bps,max_inter_len,GC_content,max_ED

    plot_historgramm('E', 'MFE', plot_dir, df_neg,df_pos)
    plot_historgramm('E_hybrid', 'E_hybrid', plot_dir, df_neg,df_pos)
    plot_historgramm('no_bps', 'no_bps', plot_dir, df_neg,df_pos)
    plot_historgramm('max_inter_len', 'max_inter_len', plot_dir, df_neg,df_pos)
    plot_historgramm('GC_content', 'GC_content', plot_dir, df_neg,df_pos)
    plot_historgramm('max_ED', 'max_ED', plot_dir, df_neg,df_pos)

    #df_pos_sub = df_pos['E']
    #df_pos_sub['data'] = 'pos'

    #df_neg_sub = df_neg['E']
    #df_neg_sub['data'] = 'neg'

    #df_data = pd.concat([df_pos_sub, df_neg_sub])

    # both historgram




if __name__ == '__main__':
    main()