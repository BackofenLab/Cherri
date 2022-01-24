#!/usr/bin/env python
import pandas as pd
import argparse
import os
import re
import time
import sys
import rrieval.lib as rl

# from peakhood import hoodlib
# from peakhood.hoodlib import extract_multicore_wrapper

__version__ = "0.1"

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~~~~~~~~~~~~~~
Check out available modes
~~~~~~~~~~~~~~~~~~~~~~~~~~

peakhood -h


~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

cd peakhood
python3 -m doctest hoodlib.py


"""


################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    # Tool description text.
    help_description = """

    Classification predicted of RNA RNA Interactions (RRI)

    """

    # Define argument parser.
    p = argparse.ArgumentParser(#add_help=False,
                                prog="cherri",
                                description=help_description)

    # Tool version.
    p.add_argument("-v", "--version", action="version",
                   version="cherri v" + __version__)

    # Add subparsers.
    subparsers = p.add_subparsers(help='Program modes')

    """
    Context extraction mode.
    """
    p_ex = subparsers.add_parser('eval',
                                  help='evaluate or classify RRI')
    p_ex.set_defaults(which='eval')
    # Add required arguments group.
    p_exm = p_ex.add_argument_group("required arguments")
    # Required arguments for evaluate.
    p_exm.add_argument("-i1", "--RRIs_table",
                       dest="RRIs_table",
                       required = True,
                       help= "table containg all rris that should be evalutated in the corrct fromat")
    p_exm.add_argument("-g", "--genome_file",
                       dest="genome_file",
                       action="store",
                       required=True,
                       help= "Genomic sequences .2bit file")
    p_exm.add_argument("-o", "--out_path",
                       required=True,
                       help= "path to folder all output folder of each step of the data preparation")
    p_exm.add_argument("-l", "--chrom_len_file",
                       action="store",
                       dest="chrom_len_file",
                       required=True,
                       help= "tabular file containing chrom name \t chrom lenght for each chromosome")


    # Optional arguments for evaluate.
    p_ex.add_argument("-i2", "--occupyed_regions",
                      default="none",
                      help= "path to occupyed regions file or specify: human, mouse or none")
    p_ex.add_argument("-c", "--context",
                      nargs='?',
                      type=int,
                      dest="context",
                      default=50,
                      help= "how much context should be added at left an right of the sequence")
    p_ex.add_argument("-n", "--experiment_name",
                      action="store",
                      dest="experiment_name",
                      default='eval_rri',
                      help= "name of the datasoruce of RRIs")
    p_ex.add_argument("-p", "--param_file",
                       dest="param_file",
                       default="../IntaRNA_param.txt",
                       help= "IntaRNA parameter file")


    """
    Build a model form new data: train
    """
    p_mrg = subparsers.add_parser('train',
                                  help='Build a model form new data: training')
    p_mrg.set_defaults(which='train')
    # Add required arguments group.
    p_mrgm = p_mrg.add_argument_group("required arguments")
    # Required arguments for merge.
    p_mrgm.add_argument("-i1", "--RRI_path",
                        required=True,
                        help= "path to folder storing the ChiRA interaction summary files",
                        default="/vol/scratch/data/RRIs/Paris/")
    #p_mrgm.add_argument("-i2", "--rbp_path",
                        #help= "path to RBP side data file (bed format)",
                        #default="/vol/scratch/data/human_RBP_coverage/GSE38355_ProtOccProf_4SU_consensus_TC_hg38.bed")
    p_mrgm.add_argument("-o", "--out_path",
                        required=True,
                        help= "path to folder all output folder of each step of the data preparation")
    p_mrgm.add_argument("-r", "--list_of_replicats",
                        action="store",
                        required=True,
                        nargs='+',
                        dest="list_of_replicats",
                        help= "list ChiRA interaction summary files names of all replicats")
    p_mrgm.add_argument("-l", "--chrom_len_file",
                        action="store",
                        dest="chrom_len_file",
                        required=True,
                        help= "tabular file containing chrom name \t chrom lenght for each chromosome")
    p_mrgm.add_argument("-g", "--genome_file",
                        action="store",
                        dest="genome_file",
                        required=True,
                        help= "path to 2bit genome file")

    # Optional arguments for merge.
    p_mrg.add_argument("-c", "--context",
                       nargs='?',
                       type=int,
                       dest="context",
                       default=50,
                       help= "how much context should be added at left an right of the sequence")
    p_mrg.add_argument("-n", "--experiment_name",
                       action="store",
                       dest="experiment_name",
                       default='model_rri',
                       help= "name of the datasoruce of RRIs")
    p_mrg.add_argument("-p", "--param_file",
                       dest="param_file",
                       help= "IntaRNA parameter file",
                       default="../IntaRNA_param.txt")

    return p


################################################################################


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




def main_eval(args):
    """

    Useful output:

    """

    print("Running for you in EXTRACT mode ... ")


    """
    Output files.

    """

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
                    ' -p ' + param_file)

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


################################################################################

def main_train(args):
    """
    Generat a model

    """
    args = parser.parse_args()
    input_path_RRIs = args.RRI_path
    #file_rbp_pos = args.rbp_path
    genome_file = args.genome_file
    replicats = args.list_of_replicats
    out_path = args.out_path
    context = args.context
    experiment_name = args.experiment_name
    chrom_len_file = args.chrom_len_file
    param_file = args.param_file

    overlap_th = 0.3

    # define output folder
    timestr = time.strftime("%Y%m%d")
    out_path = '/' + out_path + '/' + timestr + '_Cherri_model_build/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print('***added new folder***')


    ### Call find_occupied_regions.py ##########################################
    occupied_regions_param = (' -i1 ' + input_path_RRIs + ' -i2 non' +
                             ' -r ' + ' '.join(replicats) +
                             ' -o ' + out_path + ' -t ' + str(overlap_th) +
                             ' -s 0.5 ')

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
                    str(context) + ' --pos_occ -b 40  -l ' + chrom_len_file +
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


    # Convert function from biofilm to prepare data?



################################################################################

if __name__ == '__main__':
    # Setup argparse.
    parser = setup_argument_parser()
    # Print help if no parameter is set.
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
    # Read in command line arguments.
    args = parser.parse_args()


    # Are my tools ready?
    #assert hoodlib.is_tool("bedtools"), "bedtools not in PATH"
    #assert hoodlib.is_tool("twoBitToFa"), "twoBitToFa not in PATH"
    #assert hoodlib.is_tool("twoBitInfo"), "twoBitInfo not in PATH"

    # Run selected mode.
    if args.which == 'eval':
        main_eval(args)
    elif args.which == 'train':
        main_train(args)



################################################################################
