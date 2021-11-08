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

def plot_historgramm_2(colum, file_name, plot_dir, df_neg,df_pos, name, flag_norm):
    flag_flip_legend = False
    fig = plt.figure()
    if flag_norm == 1:
        n, bins, edges = plt.hist([df_neg[colum], df_pos[colum]],
                                  color=['r','b'], bins=15, alpha=0.5,
                                  density=True)
    elif flag_norm == 0:
        n, bins, edges = plt.hist([df_neg[colum], df_pos[colum]],
                                  color=['r','b'], bins=15, alpha=0.5,
                                  density=False)


    plt.xticks(bins)
    plt.xticks(rotation=75)
    #print(min(df_neg[colum]))
    #print(max([-3,-5]))
    max_val = max([max(df_neg[colum]),max(df_pos[colum])])
    min_val = min([min(df_neg[colum]),min(df_pos[colum])])
    if flag_flip_legend and min_val<0 and max_val <0:
        max_val_temp = max_val
        max_val = min_val
        min_val = max_val_temp
    plt.xlim(min_val, max_val)
    labels= ["negative", "positive"]
    plt.legend(labels)
    plt.xlabel(colum)
    plt.title(name)
    fig.savefig(plot_dir + file_name + '_histogramm.pdf', format='pdf', dpi=300, bbox_inches='tight')



def plot_historgramm_3(colum, file_name, plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm):
    flag_flip_legend = False
    fig = plt.figure()
    if flag_norm == 1:
        n, bins, edges = plt.hist([df_neg[colum], df_pos1[colum],df_pos2[colum]],
                                  color=['r','b','g'], bins=15, alpha=0.5,
                                  density=True)
    elif flag_norm == 0:
        n, bins, edges = plt.hist([df_neg[colum], df_pos1[colum],df_pos2[colum]],
                                  color=['r','b','g'], bins=15, alpha=0.5,
                                  density=False)


    plt.xticks(bins)
    plt.xticks(rotation=75)
    #print(min(df_neg[colum]))
    #print(max([-3,-5]))
    max_val = max([max(df_neg[colum]),max(df_pos1[colum]),max(df_pos2[colum])])
    min_val = min([min(df_neg[colum]),min(df_pos1[colum]),min(df_pos2[colum])])
    if flag_flip_legend and min_val<0 and max_val <0:
        max_val_temp = max_val
        max_val = min_val
        min_val = max_val_temp
    plt.xlim(min_val, max_val)
    labels= ["negative seedMaxE0", "positive seedMaxE0", "positive seedMaxE3"]
    plt.legend(labels)
    plt.xlabel(colum)
    plt.title(name)
    fig.savefig(plot_dir + file_name + '_histogramm.pdf', format='pdf', dpi=300, bbox_inches='tight')


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--pos1_file",
                        help= "file to features of positive data",
                        default="/vol/scratch/data/features_files/test/fearutes_pos.csv")
    parser.add_argument("-i2", "--neg_file",
                        help= "file to features of negative data",
                        default="/vol/scratch/data/features_files/test/fearutes_neg.csv")
    parser.add_argument("-i3", "--pos2_file",
                        help= "file to features of positive data seedMaxE3",
                        default="")
    parser.add_argument("-n", "--name",
                        help= "name to add to the output file name",
                        default='')
    parser.add_argument("-f", "--flag_norm",
                        help= "1: normalize each bin, so that the area under the histogram integrates to 1",
                        default=0)
    parser.add_argument("-o", "--out_dir",
                        help= "output directory",
                        default="/vol/scratch/data/features_files/test/")




    args = parser.parse_args()

    pos1_file = args.pos1_file
    pos2_file = args.pos2_file
    neg_file = args.neg_file
    name = args.name
    out_dir = args.out_dir
    flag_norm = int(args.flag_norm)
    # flag_norm = True
    if flag_norm == 1:
        print('histrogram will be normalized')
        file_add_name = name + '_norm_bins_'
    elif flag_norm == 0:
        print('data not normalized')
        file_add_name = name
    else:
        print('error: flag_norm should be 1 for True and 0 for False')

    plot_dir = out_dir + file_add_name + '_'


    #print(df_pos.info())
    #print(df_neg.info())

    # print E, E_hybrid, ED1, ED2, no_bps,max_inter_len,GC_content,max_ED

    if not pos2_file:
        df_pos1 = pd.read_table(pos1_file, sep=',')
        df_neg = pd.read_table(neg_file, sep=',')

        plot_historgramm_2('E', 'MFE', plot_dir, df_neg,df_pos1, name, flag_norm)
        plot_historgramm_2('E_hybrid', 'E_hybrid', plot_dir, df_neg,df_pos1, name, flag_norm)
        plot_historgramm_2('no_bps', 'no_bps', plot_dir, df_neg,df_pos1, name, flag_norm)
        plot_historgramm_2('max_inter_len', 'max_inter_len', plot_dir, df_neg,df_pos1, name, flag_norm)
        plot_historgramm_2('GC_content', 'GC_content', plot_dir, df_neg,df_pos1, name, flag_norm)
        plot_historgramm_2('max_ED', 'max_ED', plot_dir, df_neg, df_pos1, name, flag_norm)
        plot_historgramm_2('sum_ED', 'sum_ED', plot_dir, df_neg, df_pos1, name, flag_norm)
        #plot_historgramm_2('no_seeds', 'no_seeds', plot_dir, df_neg, df_pos1, name, flag_norm)
        #plot_historgramm_2('max_seed_E', 'max_seed_E', plot_dir, df_neg, df_pos1, name, flag_norm)
    else:
        df_pos1 = pd.read_table(pos1_file, sep=',')
        df_pos2 = pd.read_table(pos2_file, sep=',')
        df_neg = pd.read_table(neg_file, sep=',')

        plot_historgramm_3('E', 'MFE', plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('E_hybrid', 'E_hybrid', plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('no_bps', 'no_bps', plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('max_inter_len', 'max_inter_len', plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('GC_content', 'GC_content', plot_dir, df_neg,df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('max_ED', 'max_ED', plot_dir, df_neg, df_pos1, df_pos2, name, flag_norm)
        plot_historgramm_3('sum_ED', 'sum_ED', plot_dir, df_neg, df_pos1, df_pos2, name, flag_norm)

        # seedE,no_seeds




    #df_pos_sub = df_pos['E']
    #df_pos_sub['data'] = 'pos'

    #df_neg_sub = df_neg['E']
    #df_neg_sub['data'] = 'neg'

    #df_data = pd.concat([df_pos_sub, df_neg_sub])

    # both historgram




if __name__ == '__main__':
    main()
