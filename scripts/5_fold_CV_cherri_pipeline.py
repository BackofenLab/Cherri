#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import re
import time
import sys
import rrieval.lib as rl
import random
import json


def concat_coulums(df):
    df['ID'] = df['target_ID'] + '_' + df['query_ID']
    return df

def get_list(df):
    ID_list = df['ID'].to_list()
    unique_sorted_list = sorted(list(set(ID_list)))
    #print(len(unique_sorted_list))
    return unique_sorted_list

def split_list(in_list):
    shuffled_list = random.sample(in_list, len(in_list))
    n = len(shuffled_list) // 5
    sets = [shuffled_list[i:i+n] for i in range(0, len(shuffled_list), n)]
    return sets

def write_json_index(list_dfs, file_json):
    data_df = pd.concat(list_dfs, ignore_index=True)
    #print(data_df)
    #print(len(data_df))
    data_ID_dict = data_df['ID'].to_dict()

    print(len(data_ID_dict))
    with open(file_json, 'w') as file:
        json.dump(data_ID_dict, file)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_file",
                        help= "inputfile split into 5 sections",
                        default="/vol/scratch/data_storage/Cross_Validation/data/human/PARIS_human_context_150_pos_occ__block_ends_40_RRI_dataset.csv")
    parser.add_argument("-oc", "--occupied_regions",
                        help= "occupied regions of the added file [Format: BED]",
                        default="/vol/scratch/data/RRIs/Paris/")
    parser.add_argument("-ca", "--context_additional",
                        help= "context to extend left and right for the BED file instances",
                        default="5")
    parser.add_argument("-es", "--exp_score_th",
                        help= "score threshold for the additional occupied regions [BED]",
                        default="10")

    args = parser.parse_args()
    RRI_file = args.input_file
    #pos_file = '/vol/scratch/data_storage/Cross_Validation/data/human/PARIS_human_context_150_pos_occ_pos.csv'
    #neg_file = '/vol/scratch/data_storage/Cross_Validation/data/human/PARIS_human_context_150_pos_occ_neg.csv'

    pos_file = '/vol/scratch/data_storage/Cross_Validation/data/human/test2_context_50_pos_occ_pos.csv'
    neg_file = '/vol/scratch/data_storage/Cross_Validation/data/human/test2_context_50_pos_occ_neg.csv'
    file_json = '/vol/scratch/data_storage/Cross_Validation/data/human/index_dict.json'

    pos_df = rl.read_chira_data(pos_file, 'yes', ',')
    pos_df = concat_coulums(pos_df)

    neg_df = rl.read_chira_data(neg_file, 'yes', ',')
    neg_df = concat_coulums(neg_df)

    print(len(pos_df))
    print(len(neg_df))

    write_json_index([pos_df,neg_df], file_json)
    index_values = neg_df.index.tolist()
    print(index_values)




    pos_IDs_list = pos_df['ID'].to_list()
    neg_IDs_list = neg_df['ID'].to_list()

    pos_unique_sorted_list = get_list(pos_df)
    neg_unique_sorted_list = get_list(neg_df)

    overlap = list(set(pos_unique_sorted_list).intersection(neg_unique_sorted_list))
    #print(len(overlap))
    #print(overlap)

    list_cv_split = split_list(overlap)
    #print(list_cv_split)


    # Daten nach liste sotiren










# read data: from the occypied regions

# solve issue of th! the bed files have a score_th

#bed files aus den trusted RRIs machen
#wenn RBP mit RBP concatinieren
#Cherri starten
#Outputs verwalten
#Evaluation starten
#TP’s und FN’s berechnen
#Ergebniss speichern
#Über alle 5 files loopen
#Ergebniss zusammen fügen




if __name__ == '__main__':
    main()
