#!/usr/bin/env python
import pandas as pd
import argparse
import os
import time
import rrieval.lib as rl



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--RRI_path",
                        help= "path to folder storing the ChiRA interaction summary files",
                        default="/vol/scratch/data/RRIs/Paris/")
    parser.add_argument("-i2", "--rbp_path",
                        help= "path to RBP side data file (bed format)",
                        default="/vol/scratch/data/human_RBP_coverage/GSE38355_ProtOccProf_4SU_consensus_TC_hg38.bed")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file",
                        required=True, help= "path to 2bit genome file")
    parser.add_argument("-r", "--list_of_replicats", action="store",
                        nargs='+',
                        dest="list_of_replicats", required=True,
                        help= "list ChiRA interaction summary files names of all replicats")
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


    args = parser.parse_args()
    input_path_RRIs = args.RRI_path
    file_rbp_pos = args.rbp_path
    genome_file = args.genome_file
    replicats = args.list_of_replicats
    out_path = args.out_path
    context = args.context
    experiment_name = args.experiment_name
    chrom_len_file = args.chrom_len_file

    overlap_th = 0.3

    # define output folder
    timestr = time.strftime("%Y%m%d")
    out_path = '/' + out_path + '/' + timestr + '_Cherri_model_build/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print('***added new folder***')


    ### Call find_occupied_regions.py ##########################################
    occupied_regions_param = (' -i1 ' + input_path_RRIs + ' -i2 ' +
                             file_rbp_pos + ' -r ' + ' '.join(replicats) +
                             ' -o ' + out_path + ' -t ' + str(overlap_th))

    call_occupied_regions = ('python -W ignore find_occupied_regions.py' +
                            occupied_regions_param)
    print('occupyed regions call:')
    print(call_occupied_regions)
    rl.call_script(call_occupied_regions)

    # output:
    occupyed_outfile =  out_path + '/' + timestr + '_occ_out/'



    #### Call find_occupied_regions.py #########################################
    pos_neg_out_path = out_path + 'pos_neg_data/'
    os.mkdir(pos_neg_out_path)
    trusted_rri_file = (occupyed_outfile + 'rri_occupied_regions_overlapTH_' +
                       str(overlap_th) + '_scoreTH_1.cvs')
    occupyed_regions_file =  occupyed_outfile + 'occupied_regions.obj'

    pos_neg_param = (' -i1 ' + trusted_rri_file + ' -i2 ' +
                    occupyed_regions_file + ' -d ' + pos_neg_out_path + ' -g ' +
                    genome_file + ' -n ' + experiment_name + ' -c ' +
                    str(context) + ' --pos_occ -b 40  -l ' + chrom_len_file)

    call_pos_neg = ('python -W ignore generate_pos_neg_with_context.py' +
                            pos_neg_param)
    print('Positive and negative instaces call:')
    print(call_pos_neg)
    rl.call_script(call_pos_neg)


    ### Call get_features.py ###################################################
    feature_out_path = out_path + 'feature_files/'
    os.mkdir(feature_out_path)
    midel_name =  experiment_name + '_context_' + str(context)
    file_neg = pos_neg_out_path + midel_name + '_pos_occ_neg.csv'
    file_pos = pos_neg_out_path + midel_name + '_pos_occ_pos.csv'
    feature_pos = (feature_out_path +  'feature_filtered_' + midel_name +
                  '_pos_occ_pos.csv')
    feature_neg = (feature_out_path +  'feature_filtered_' + midel_name +
                  '_pos_occ_neg.csv')

    pos_feature_param = ' -i ' + file_pos + ' -f all -o ' + feature_pos
    call_pos_feature = 'python -W ignore get_features.py' + pos_feature_param
    print('Positive feature call:')
    print(call_pos_feature)
    rl.call_script(call_pos_feature)

    neg_feature_param = ' -i ' + file_neg + ' -f all -o ' + feature_neg
    call_neg_feature = 'python -W ignore get_features.py' + neg_feature_param
    print('Negative feature call:')
    print(call_neg_feature)
    rl.call_script(call_neg_feature)


if __name__ == '__main__':
    main()
