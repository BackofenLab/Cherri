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
    parser.add_argument("-s", "--input_set", action="store", dest="input_set", required=True
                        , help= "prefix indicating the experiment and subest to be analyzed")
    parser.add_argument("-f", "--feature_set_list", action="store", nargs='+',
                        dest="feature_set_list", required=True,
                        help= "set of features the script will output ")
    parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True
                        , help= "dir where feature tables will be stored")




    args = parser.parse_args()

    input_dir = args.input_dir
    input_set = args.input_set
    feature_set_list = args.feature_set_list
    output_dir = args.output_dir

    print(feature_set_list)
    #output_dir = '/vol/scratch/data/feature_input_HEK293T/'


    calls_dict = {}
    overview_dict = {}
    counter = 1
    #file_names_list = ['paris_mES_06_context_method_together_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_3_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_2_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_3_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_2_with_10_context_']

    file_names_list_HEK293T = ['paris_HEK293T_context_method_together_shuffling_method_3_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_2_with_10_context_', 'paris_HEK293T_context_method_together_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_1_with_10_context_', 'paris_HEK293T_context_method_separat_shuffling_method_2_with_10_context_']
    file_names_list_mES = ['paris_mES_06_context_method_together_shuffling_method_3_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_2_with_10_context_', 'paris_mES_06_context_method_together_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_1_with_10_context_', 'paris_mES_06_context_method_separat_shuffling_method_2_with_10_context_']
    file_names_list_lncRNA = ['paris_lncRNA_context_method_together_shuffling_method_3_with_10_context_', 'paris_lncRNA_context_method_together_shuffling_method_2_with_10_context_', 'paris_lncRNA_context_method_together_shuffling_method_1_with_10_context_', 'paris_lncRNA_context_method_separat_shuffling_method_1_with_10_context_', 'paris_lncRNA_context_method_separat_shuffling_method_2_with_10_context_']
    file_names_list_snRNA = ['paris_snRNA_context_method_together_shuffling_method_3_with_10_context_', 'paris_snRNA_context_method_together_shuffling_method_2_with_10_context_', 'paris_snRNA_context_method_together_shuffling_method_1_with_10_context_', 'paris_snRNA_context_method_separat_shuffling_method_1_with_10_context_', 'paris_snRNA_context_method_separat_shuffling_method_2_with_10_context_']

    if input_set == 'human':
        file_names_list = file_names_list_HEK293T
    elif input_set == 'mouse':
        file_names_list = file_names_list_mES
    elif input_set == 'lncRNA':
        file_names_list = file_names_list_lncRNA
    elif input_set == 'snRNA':
        file_names_list = file_names_list_snRNA



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
    f1.write('#!/usr/bin/env bash\n\n')
    #f1.write('\n\ntrap ctrl_c INT\n\n\nfunction ctrl_c() {\necho "** Trapped CTRL-C"\nexit\n}\n\n\n\n')
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
