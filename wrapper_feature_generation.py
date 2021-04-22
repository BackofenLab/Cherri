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
    output_dir = '/vol/scratch/data/feature_input_HEK293T/'


    calls_dict = {}
    overview_dict = {}
    counter = 1
    #file_names_list = ['paris_mES_06_context_method_together_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_3_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_2_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_3_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_2_with_10_context_']

    file_names_list = ['paris_HEK293T_context_method_together_shuffling_method_3_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_2_with_10_context_']

    print(len(file_names_list))


    for file_name in file_names_list:
        pos_file = input_dir + file_name + '_pos_RRI_dataset.csv'
        neg_file = input_dir + file_name + '_neg_RRI_dataset.csv'
        for featurs in feature_set_list:
            id_pos = str(counter)+'a'
            id_neg = str(counter)+'b'
            call_pos = 'python get_features.py -i ' + pos_file + ' -f ' + featurs + ' -o ' + output_dir + id_pos
            call_neg = 'python get_features.py -i ' + neg_file + ' -f ' + featurs + ' -o ' + output_dir + id_neg
            overview_pos = id_pos + ' \t ' + file_name + ' \t ' + ' pos \t ' + featurs + '\n'
            overview_neg = id_neg + ' \t ' + file_name + ' \t ' + ' neg \t ' + featurs + '\n'
            calls_dict[id_pos]= call_pos
            calls_dict[id_neg]= call_neg
            overview_dict[id_pos]= overview_pos
            overview_dict[id_neg]= overview_neg
            counter+= 1
    #print(overview_dict)

    f1 = open(output_dir + "calls.sh", "w")
    f2 = open(output_dir + "overview.tabular", "w")
    # write bash script!
    f1.write('#!/usr/bin/env bash\n\ntrap ctrl_c INT\n\n\n')
    f1.write('function ctrl_c() {\necho "** Trapped CTRL-C"\nexit\n}\n\n\n\n')
    # header table
    f2.write('id\tsuffeling_id\tdata\tfeaturs\n')
    for key in calls_dict:
        f1.write('# call %s\n'%(key))
        f1.write(calls_dict[key] + '\n')
        f2.write(overview_dict[key])
    f1.close()
    f2.close()



if __name__ == '__main__':
    main()
