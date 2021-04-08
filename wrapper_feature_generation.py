#!/usr/bin/env python
import pandas as pd
import math
import matplotlib as mpl
import argparse
import rrieval.lib as rl





def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True
                                           , help= "path to input files")
    parser.add_argument("-f", "--feature_set_list", action="store", nargs='+',
                        dest="feature_set_list", required=True,
                        help= "set of features the script will output ")



    args = parser.parse_args()

    input_dir = args.input_dir
    feature_set_list = args.feature_set_list

    print(feature_set_list)
    output_dir = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/input_features_test/'


    calls_dict = {}
    overview_dict = {}
    counter = 1
    file_names_list = ['paris_mES_06_context_method_together_shuffling_method_2_with_10_context', 'paris_HEK293T_06_context_method_together_shuffling_method_2_with_10_context']

    for file_name in file_names_list:
        pos_file = input_dir + file_name + '_pos_RRI_dataset.csv'
        neg_file = input_dir + file_name + '_neg_RRI_dataset.csv'
        for featurs in feature_set_list:
            id_pos = str(counter)+'a'
            id_neg = str(counter)+'b'
            call_pos = 'python get_features.py -i ' + pos_file + ' -f ' + featurs + ' -o ' + output_dir + id_pos + '/'
            call_neg = 'python get_features.py -i ' + neg_file + ' -f ' + featurs + ' -o ' + output_dir + id_neg + '/'
            overview_pos = id_pos + ' \t ' + file_name + ' \t ' + ' pos \t ' + featurs + '\n'
            overview_neg = id_neg + ' \t ' + file_name + ' \t ' + ' neg \t ' + featurs + '\n'
            calls_dict[id_pos]= call_pos
            calls_dict[id_neg]= call_neg
            overview_dict[id_pos]= overview_pos
            overview_dict[id_neg]= overview_pos
            counter+= 1
    print(overview_dict)

    f1 = open(output_dir + "calls.txt", "w")
    f2 = open(output_dir + "overview.tabular", "w")
    f2.write('id  \t suffeling_id  \t data  \t featurs\n')
    for key in calls_dict:
        f1.write(calls_dict[key] + '\n')
        f2.write(overview_dict[key])
    f1.close()
    f2.close()



if __name__ == '__main__':
    main()
