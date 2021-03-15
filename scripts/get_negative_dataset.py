#!/usr/bin/env python
import pandas as pd
import math
#import matplotlib as mpl
from collections import defaultdict
import sys
#import seaborn
import argparse
#import csv
#import collections
#from operator import itemgetter
#import sys
import subprocess
import random
import re
import os
import time


def shuffle_sequence(seq, times):
    """
    shuffle on given sequence x times

        Parameters
        ----------
        seq: sequence
        times: amount of shuffeling

        Raises
        ------
        nothing

        Returns
        -------
        seq_list
            list of shuffled sequences

        """
    #seq_list= []

    call = "ushuffle -s " + str(seq) + " -n " \
                       + str(times) + " -k 2"
    #print(call)
    p = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True,)
    stdout, stderr = p.communicate()
    seq_list = stdout.decode('utf-8').strip().split('\n')

    return seq_list

def get_neg_instance(df_pos, df_neg):
    """
    Initializing dataframe with the IntaRNA coulum header

        Parameters
        ----------
        df_pos: dataframe containing the positive RRI InraRNA prediction
        df_neg: dataframe conainging the IntrRNA reusults of the suffled RRI
        sequences of the pos RRI

        Raises
        ------
        nothing

        Returns
        -------
        df_neg_entry
            datafram containing the one neg instance with the closed enery
            to the positive instance

        """

    E_pos_RRI = float(df_pos.E.tolist()[0])
    E_neg_RRI_list = [float(i) for i in df_neg.E.tolist()]

    closest_E = min(E_neg_RRI_list, key=lambda x:abs(x-E_pos_RRI))

    df_neg_entry = df_neg[df_neg['E'] == str(closest_E)]
    #print(df_neg_entry['query'])
    return df_neg_entry



def inial_df():
    """
    Initializing dataframe with the IntaRNA coulum header

        Parameters
        ----------


        Raises
        ------
        nothing

        Returns
        -------
        df_inital
            emty df with just the header (coloum)

        """
    df_inital = pd.DataFrame(columns=['id1','start1','end1','id2','start2',
                                      'end2','subseqDP','hybridDP','E', 'target', 'query'])
    return df_inital


def Intarna_call(seq1, seq2,df):
    """
    Intarna_call

        Parameters
        ----------
        seq1: tartet sequence
        seq2: query sequence
        df: datafame where the IntaRNA result will be appended

        Raises
        ------
        nothing

        Returns
        -------
        df_result
            df incuding IntaRNA result

        """
    #print(df)
    temp_out_cvs_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/test_IntaRNA/test.csv'
    #call = 'IntaRNA -t ' + seq1 + ' -q ' + seq2 + ' --out ' + temp_out_cvs_path + ' --outMode C'
    call = 'IntaRNA -t ' + seq1 + ' -q ' + seq2 + ' --outMode C'
    print(call)
    process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
    for idx, line in enumerate(process.stdout):
        line = line.decode("utf-8").strip().split(';')
        #line = line.strip()
        #line = line.split(';')
        print(line)
        if idx == 0:
            col = line
            result = 'empty'
        elif idx == 1:
            values = line
            result = 'exists'
        else:
            print('error IntaRNA output has unexpected number of rows')
        #print(line)
    if result == 'empty':
        values = ['nan','nan','nan','nan','nan','nan','nan','nan','nan', 'nan', 'nan'])
        df_one_result = pd.DataFrame([values], columns=col)

    df_one_result = pd.DataFrame([values], columns=col)
    df_one_result['target'] = seq1
    df_one_result['query'] = seq2
    df_result = pd.concat([df, df_one_result])

    #time.sleep(15)
    #assert os.path.isfile(temp_out_cvs_path), 'temp_out_cvs_path not found'
    #if os.stat(temp_out_cvs_path).st_size == 0:
        #print('file is empty')
        #print(temp_out_cvs_path)
        #sys.exit()

    #f = open(temp_out_cvs_path, "r")
    #print(f.read())
    #df_one_result = pd.read_table(temp_out_cvs_path, sep=";")
    #df_result = pd.concat([df, df_one_result])
    #os.remove(temp_out_cvs_path)
    #else:
        #df_result = df
    return df_result



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                        required=True,
                        help= "path to file storing all positve trusted RRIs")


    args = parser.parse_args()
    input_file = args.input_file

    df_pos_RRIs_result = inial_df()
    df_neg_RRIs_result = inial_df()

    # df.columns.tolist()
    df_RRIs = pd.read_table(input_file, sep=",")
    target_seq_list = df_RRIs.ineraction_side_1st.tolist()
    query_seq_list = df_RRIs.ineraction_side_2end.tolist()
    #seq1_list = ['GAGCUCCCGGGGGGGGGGGGGGGGGGGGGGGGCCA']
    #seq2_list = ['AAAACCCCCCCUUUU']

    test_no_seq = 5
    for idx, target in enumerate(target_seq_list):
        query = query_seq_list[idx]
        # inital df for output
        df_initial_pos_result = inial_df()
        # intarna call to get pos enegy
        df_pos = Intarna_call(target, query, df_initial_pos_result)

        # temp neg data
        shuffled_target_list = shuffle_sequence(target, test_no_seq)
        #print(shuffled_target_list)
        shuffled_query_list = shuffle_sequence(query, test_no_seq)
        #print(shuffled_query_list)
        df_initial_neg_result = inial_df()
        for idx2, neg_target in enumerate(shuffled_target_list):
            neg_query = shuffled_query_list[idx2]
            if idx2 == 0:
                #print('if')
                df_neg = Intarna_call(neg_target, neg_query, df_initial_neg_result)
            else:
                #print('else')
                df_neg = Intarna_call(neg_target, neg_query, df_neg)
        #print(df_pos)
        #print(df_neg)
        df_neg_entry = get_neg_instance(df_pos, df_neg)
        df_pos_RRIs_result = pd.concat([df_pos_RRIs_result, df_pos])
        df_neg_RRIs_result = pd.concat([df_neg_RRIs_result, df_neg_entry])


if __name__ == '__main__':
    main()
