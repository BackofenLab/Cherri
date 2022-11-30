#!/usr/bin/env python3
import pandas as pd
import math
import matplotlib as mpl
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import os
#import os.path
import argparse
from matplotlib.patches import Rectangle
import rrieval.lib as rl

def include_zero_rows(df1, col1, df2, col2, hue):
    """
    check which RNAs are not in the other df and append them with a count of 0

        Parameters
        ----------
        df1: dataframe 1
        col1: colum name of the first df (RNA)
        df1: dataframe 2
        col2: colum name of the second df (RNA)

        Raises
        ------
        nothing

        Returns
        -------
        df1
            df incuding the added RNAs with 0 count

        """
    #print(df1[col1])
    #print(df2[col2])
    no_existing_RNAs = list(set(df1[col1])-set(df2[col2]))
    #print(no_existing_RNAs)
    for rna in no_existing_RNAs:
        row= [rna, 0, hue]
        df2.loc[len(df2), :] = row
    return df1


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--input_pos", action="store", dest="input_pos", required=True
                                           , help= "positive data table file.")
    parser.add_argument("-i2", "--input_neg", action="store", dest="input_neg", required=True
                                           , help= "negative data table file")
    parser.add_argument("-i3", "--input_dir", action="store", dest="input_dir", required=True
                                           , help= "path to the input data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True
                                      , help= "path to save data to.")

    args = parser.parse_args()

    input_pos = args.input_pos
    input_neg = args.input_neg
    input_dir = args.input_dir
    plot_dir = args.save_path

    pos_file = input_dir + input_pos
    neg_file = input_dir + input_neg

    df_pos = rl.read_chira_data(pos_file, header='yes', separater=",")
    df_neg = rl.read_chira_data(neg_file, header='yes', separater=",")
    # print(df_pos.info())

    plot_dir = plot_dir + 'plots/'
    if os.path.isdir(plot_dir):
        print('plot folder already exists: %s'%(plot_dir))
    else:
        os.mkdir(plot_dir)

    ####RNAs#######################################################
    df_RNA_first = df_neg['biotype_region_1st'].value_counts().reset_index()
    df_RNA_first.columns = ['RNA', 'count']
    df_RNA_first['hue'] ='first'

    df_RNA_second = df_neg['biotype_region_2end'].value_counts().reset_index()
    df_RNA_second.columns = ['RNA', 'count']
    df_RNA_second['hue'] ='second'

    df_first =  include_zero_rows(df_RNA_first, 'RNA', df_RNA_second, 'RNA', 'second')
    df_second =  include_zero_rows(df_RNA_second, 'RNA', df_RNA_first, 'RNA', 'first')

    df_RNA = pd.concat([df_first, df_second])

    sns.set(style="whitegrid", font_scale=1)
    fig, ax = plt.subplots()
    sns.barplot(x="RNA", y="count", data=df_RNA, hue="hue", ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.savefig(plot_dir + 'RNA_his_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')


    ####Energy##########################################################

    # neg data
    sns.set(style="whitegrid", font_scale=1)
    fig, ax = plt.subplots()
    sns.distplot(df_neg['E'], ax=ax)
    fig.savefig(plot_dir + 'E_neg_distplot.pdf', format='pdf', dpi=300, bbox_inches='tight')

    # pos data
    sns.set(style="whitegrid", font_scale=1)
    fig, ax = plt.subplots()
    sns.distplot(df_pos['E'])
    fig.savefig(plot_dir + 'E_pos_distplot.pdf', format='pdf', dpi=300, bbox_inches='tight')

    # both historgram
    fig = plt.figure()
    n, bins, edges = plt.hist([df_neg['E'], df_pos['E']], color=['r','b'], bins=10, alpha=0.5)
    plt.xticks(bins)
    plt.xticks(rotation=75)
    labels= ["negative", "positive"]
    plt.legend(labels)
    fig.savefig(plot_dir + 'E_his_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')


    ####interaction length######################################################

    df_neg['avg_len_interaction'] = ((df_neg['end1'] - df_neg['start1']) + (df_neg['end2'] - df_neg['start2']))/2
    df_pos['avg_len_interaction'] = ((df_pos['end1'] - df_pos['start1']) + (df_pos['end2'] - df_pos['start2']))/2

    fig = plt.figure()
    n, bins, edges = plt.hist([df_neg['avg_len_interaction'], df_pos['avg_len_interaction']], color=['r','b'], bins=10, alpha=0.5)
    plt.xticks(bins)
    plt.xticks(rotation=75)
    labels= ["negative", "positive"]
    plt.legend(labels)
    fig.savefig(plot_dir + 'len_interaction_his_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    main()
