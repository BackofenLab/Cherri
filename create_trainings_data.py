#!/usr/bin/env python
import pandas as pd
from collections import defaultdict
import argparse
from interlap import InterLap
import subprocess
import os
import time
import rrieval.lib as rl
import pickle



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--RRI_path",
                        help= "path to folder storing all RRI data (tabel)",
                        default="/vol/scratch/data/RRIs/Paris/")
    parser.add_argument("-i2", "--rbp_path",
                        help= "path to RBP side data file (bed format)",
                        default="/vol/scratch/data/human_RBP_coverage/GSE38355_ProtOccProf_4SU_consensus_TC_hg38.bed")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file",
                        required=True, help= "path to 2bit genome file")
    parser.add_argument("-r", "--list_of_replicats", action="store",
                        nargs='+',
                        dest="list_of_replicats", required=True,
                        help= "list having filenames of all replicats")
    parser.add_argument("-o", "--out_path",
                        help= "path to folder storing outputfiles",
                        default="/vol/scratch/data/RRIs/")
    parser.add_argument("-c", "--context",  nargs='?', type=int,
                        dest="context",  default=5,
                        help= "how much context should be added at left an right of the sequence")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of positve trusted RRIs")
    parser.add_argument("-l", "--chrom_len_file",  action="store", dest="chrom_len_file",
                        required=True,
                        help= "tabular file containing chrom name \t chrom lenght for each chromosome")



    args = parser.parse_args()
    input_path_RRIs = args.RRI_path
    file_rbp_pos = args.rbp_path
    genome_file = args.genome_file
    replicats = args.list_of_replicats
    out_path = args.out_path
    context = args.context
    experiment_name = args.experiment_name
    chrom_len_file = args.chrom_len_file

    timestr = time.strftime("%Y%m%d")
    out_path = '/' + out_path + '/' + timestr + 'Cherri_model_build/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print('***added new folder***')


    # Call find_occupied_regions.py
    occupied_regions_param = (' -i1 ' + input_path_RRIs + ' -i2 ' +
                             file_rbp_pos + ' -r ' + replicats + ' -o ' +
                             out_path)

    call_occupied_regions = ('python -W ignore find_occupied_regions.py' +
                            occupied_regions_param)

    rl.call_script(call_occupied_regions)

    # output:
    occupyed_outfile =  out_path + '/' + timestr + '_occ_out/'



    # Call find_occupied_regions.py
    pos_neg_out_path = out_path + 'pos_neg_data/'

    pos_neg_param = (' -i1 ' + trusted_rri_file + ' -i2 ' +
                             occupyed_regions_file + ' -d ' + out_path +
                             ' -g ' + genome_file + ' -n ' + experiment_name +
                             ' -c ' + context + ' --pos_occ -b 40  -l ' +
                             chrom_len_file)

    call_pos_neg = ('python -W generate_pos_neg_with_context.py' +
                            pos_neg_param)

    rl.call_script(call_pos_neg)


    # Call get_features.py
    feature_out_path = out_path + 'feature_files/'
    midel_name =  experiment_name + '_context_' + context
    file_neg = pos_neg_out_path + midel_name + '_pos_occ_neg.csv'
    file_pos = pos_neg_out_path + midel_name + '_pos_occ_pos.csv'
    feature_pos = feature_out_path +  'feature_filtered_' + midel_name + '_pos_occ_pos.csv'
    feature_neg = feature_out_path +  'feature_filtered_' + midel_name + '_pos_occ_neg.csv'

python get_features.py
-i /vol/scratch/data/pos_neg_data_context/mouse_full_c150/paris_mouse_context_150_pos_occ_neg.csv
-f all
-o /vol/scratch/data/features_files/full_mouse_c150_maxloop3/feature_filtered_paris_mouse_context_150_pos_occ_neg.csv


    pos_feature_param = (' -i ' + file_pos + ' -f all -o ' + feature_pos)
    call_pos_feature = ('python -W get_features.py' + pos_feature_param)

    rl.call_script(call_pos_feature)

    neg_feature_param = (' -i ' + file_neg + ' -f all -o ' + feature_neg)
    call_neg_feature = ('python -W get_features.py' + neg_feature_param)

    rl.call_script(call_neg_feature)





if __name__ == '__main__':
    main()
