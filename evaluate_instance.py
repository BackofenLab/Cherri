#!/usr/bin/env python
import pandas as pd
import argparse
import os
import re
import time
import rrieval.lib as rl

def read_RRI_table(file):
    """
    Read in Table with rri interaction having the data in the format:
    chrom1,start1,stop1,stand1,chrom2,start2,stop2,strand2
        Parameters
        ----------
        file: file location of the table file

        Returns
        -------
        RRI_dict:
            interaction_id -> [chrom1,start1,stop1,stand1,chrom2,start2,stop2,strand2]

    """
    RRI_dict = {}

    # Open FASTA either as .gz or as text file.
    if re.search(".+\.gz$", file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file, "r")

    for idx, line in enumerate(f):
        if idx == 0:
            continue
        positon_list = line. rstrip("\n").split(",")
        RRI_dict[idx] = positon_list
    f.close()

    return RRI_dict


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--RRIs_table",
                        help= "table containg all rris that should be evalutated in the corrct fromat",
                        default="/vol/scratch/data/RRIs/test/test_evalueat_rris.cvs")
    parser.add_argument("-i2", "--occupyed_regions",
                        help= "path to occupyed regions file or specify: human, mouse or none",
                        default="/vol/scratch/data/human_RBP_coverage/GSE38355_ProtOccProf_4SU_consensus_TC_hg38.bed")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file",
                        required=True, help= "path to 2bit genome file")
    parser.add_argument("-o", "--out_path",
                        help= "path to folder all output folder of each step of the data preparation",
                        default="/vol/scratch/data/RRIs/")
    parser.add_argument("-c", "--context",  nargs='?', type=int,
                        dest="context",  default=5,
                        help= "how much context should be added at left an right of the sequence")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of RRIs")
    parser.add_argument("-l", "--chrom_len_file",  action="store", dest="chrom_len_file",
                        required=True,
                        help= "tabular file containing chrom name \t chrom lenght for each chromosome")
    parser.add_argument("-p", "--param_file",
                        help= "IntaRNA parameter file",
                        default="./IntaRNA_param.txt")


    args = parser.parse_args()
    RRIs_table = args.RRIs_table
    occupyed_regions = args.occupyed_regions

    genome_file = args.genome_file
    out_path = args.out_path
    context = args.context
    experiment_name = args.experiment_name
    chrom_len_file = args.chrom_len_file
    param_file = args.param_file

    overlap_th = 0.3

    # define output folder
    timestr = time.strftime("%Y%m%d")
    out_path = '/' + out_path + '/' + timestr + '_Cherri_evaluating_RRIs/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print('***added new folder***')


    ### get input data ##########################################
    # 34
    header = ['#reads','chrom_1st','start_1st','end_1st', 'strand_1st',
              'chrom_2end','start_2end','end_2end', 'strand_2end',
              'ineraction_side_1st', 'ineraction_side_2end',
              'IntaRNA_prediction', 'energy',
              'seq_1st_ineraction_side', 'seq_2end_ineraction_side',
              'start_interaction',
              'chrom_seq_1st_side', 'start_seq_1st_side',
              'stop_seq_1st_side','strand_seq_1st_side',
              'chrom_seq_2end_side', 'start_seq_2end_side',
              'stop_seq_2end_side','strand_seq_2end_side',
              'TPM_seq_1st_side', 'TPM_seq_2end_side', 'TPM_summary',
              'score_seq_1st_side', 'score_seq_2end_side','score_product',
              'biotype_region_1st', 'biotype_region_2end', 'ID_1st','ID_2end']
    end_data = ['nan', 'nan', 'nan', 'nan', 'nan','nan','nan','nan','nan','nan',
                'nan','nan', 'nan','nan','nan','nan', 'nan', 'nan','nan', 'nan',
                'nan','nan', 'nan', 'nan','nan']

    RRI_dict = read_RRI_table(RRIs_table)

    df_content_list = []
    for key in RRI_dict:
        pos_data = RRI_dict[key]
        #print(pos_data)
        # 1 + 8 + 25
        data =  ['nan'] + RRI_dict[key] + end_data
        #print(data)
        df_content_list.append(data)

    df_rris = pd.DataFrame(df_content_list,columns=header)
    eval_rri_file = out_path + 'evaluete_RRIs.table'
    df_rris.to_csv(eval_rri_file, index=False, sep=',')





    ### Call find_occupied_regions.py ##########################################
    if occupyed_regions == 'human':
        occupyed_regions = '/vol/scratch/data/RRIs/20210907_occ_out_rbp_c5/occupied_regions.obj'
    elif occupyed_regions == 'mouse':
        occupyed_regions = '/vol/scratch/data/RRIs/20210923_occ_out/occupied_regions.obj'
    elif occupyed_regions == 'none':
        occupyed_regions = 'none'
    else:
        print('using your own occupyed regions objects')


    #### Call pos data call #########################################
    pos_neg_out_path = out_path + 'positive_instance/'
    os.mkdir(pos_neg_out_path)

    pos_neg_param = (' -i1 ' + eval_rri_file + ' -i2 ' +
                    occupyed_regions + ' -d ' + pos_neg_out_path + ' -g ' +
                    genome_file + ' -n ' + experiment_name + ' -c ' +
                    str(context) + ' --no_pos_occ -s 1 -l ' + chrom_len_file +
                    ' - p ' + param_file)

    call_pos_neg = ('python -W ignore generate_pos_neg_with_context.py' +
                            pos_neg_param)
    print('Positive and negative instaces call:')
    print(call_pos_neg)
    rl.call_script(call_pos_neg)


    ### Call get_features.py ###################################################
    feature_out_path = out_path + 'feature_files/'
    os.mkdir(feature_out_path)
    midel_name =  experiment_name + '_context_' + str(context)
    file_pos = pos_neg_out_path + midel_name + 'pos.csv'
    feature_pos = (feature_out_path +  'feature_filtered_' + midel_name +
                  '_pos.csv')


    pos_feature_param = ' -i ' + file_pos + ' -f all -o ' + feature_pos
    call_pos_feature = 'python -W ignore get_features.py' + pos_feature_param
    print('Positive feature call:')
    print(call_pos_feature)
    rl.call_script(call_pos_feature)


if __name__ == '__main__':
    main()
