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



def build_interlap_occ_sides(df_interactions, flag):
    """
    Building the inerlap objects for a fast overlap comparision for one replicat

        Parameters
        ----------
        df_interactions : df including the filtered RRIs
        flag : one for rbps (only one sequences stretch in df),
               two for rris (two interacting sequences stretches in df)

        Returns
        -------
        inter_rep
            inerlap objects for a fast overlap comparision

        """

    flag_name = 'crl'
    # use defaultdict to key by chromosome.
    inter_rep = defaultdict(InterLap)

    if flag_name == 'hybrid':
        names_first = ['chrom_seq_1st_site', 'strand_seq_1st_site', 'start_seq_1st_site', 'stop_seq_1st_site']
        names_second = ['chrom_seq_2end_site', 'strand_seq_2end_site', 'start_seq_2end_site', 'stop_seq_2end_site']
    elif flag_name == 'crl':
        names_first = ['chrom_1st','strand_1st','start_1st','end_1st' ]
        names_second = ['chrom_2end','strand_2end','start_2end','end_2end']

    if flag == 'two':
        print(df_interactions['chrom_1st'])
        list_chrom_no_int = rl.get_list_chrom(df_interactions)
        df_interactions[names_first[0]] = df_interactions[names_first[0]].apply(lambda x: rl.check_convert_chr_id(x))
        df_interactions[names_second[0]] = df_interactions[names_second[0]].apply(lambda x: rl.check_convert_chr_id(x))
        print(df_interactions['chrom_1st'])
        for index, row in df_interactions.iterrows():
            row_content = row
            if not row[names_first[0]]:
                print(f'Chrom 1: {row[names_first[0]]}')
                print('can not use chromosome')
            else:
                both_keys1 = str(row[names_first[0]]) + ';' + row[names_first[1]]
                #print(f'########\n{row[names_first[3]]}\n###########')
                inter_rep[both_keys1].add((int(row[names_first[2]]), int(row[names_first[3]]), [row]))

            if not row[names_second[0]]:
                print(f'Chrom 2: {row[names_second[0]]}')
                print('can not use chromosome')
            else:
                both_keys2 = str(row[names_second[0]]) + ';' + row[names_second[1]]
                inter_rep[both_keys2].add((int(row[names_second[2]]), int(row[names_second[3]]), [row]))
    elif flag == 'one':
        list_chrom_no_int = rl.get_chrom_list_no_numbers(df_interactions, 'chrom')
        for index, row in df_interactions.iterrows():
            row_content = row

            both_keys = str(row['chrom']) + ';' + row['strand']
            inter_rep[both_keys].add((row['start'], row['end'], [row]))

    return inter_rep



def count_entrys(inter_obj, name):
    """
    count enteys of inter lap object and prints the the counts

        Parameters
        ----------
        inter_obj : inter lap object
        name: name of the data in object

    """
    count = 0
    for key in inter_obj:
        count += len(list(inter_obj[key]))
        # print(key)
        # print(list(inter_rri[key]))
    print('##########')
    print('entrys in list ',name, ' are: ' , count)


def get_prot_occ_regions(file_rbp_pos, exp_score_th, context):
    """
    get_prot_occ_regions

        Parameters
        ----------
        file_rbp_pos: bed file storing singel interaction positions
        exp_score_th: threshold to filter sequences with bed file
        context: number of nucleotide postions added to the interaction position

        Returns
        -------
        inter_rbp
            interLab object storing all cromosomal RBP sided {chrom:strand}->[s,e,rbp]


    """
    header = ['chrom', 'start', 'end', 'info', 'score', 'strand']
    df_bed_temp = pd.read_table(file_rbp_pos, header=None, sep="\t")
    df_bed = pd.DataFrame(df_bed_temp.values, columns=header)

    # filter by sorcre
    #print(df_bed)
    df_bed = df_bed[df_bed.score >= exp_score_th]
    #print(df_bed)

    # check that chorm starts with chr
    df_bed['chrom'] = df_bed['chrom'].apply(lambda x: rl.check_convert_chr_id(x))
    # add context
    df_context =  rl.add_context(df_bed, context, 'start', 'end')
    #print(df_context.info())

    inter_rep_one = build_interlap_occ_sides(df_context, 'one')
    inter_rbp = rl.mearge_overlaps(inter_rep_one, 'rbp')
    return inter_rbp


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--RRI_path",
                        help= "path to folder storing all RRI data (tabel)",
                        default="/vol/scratch/data/RRIs/Paris/")
    parser.add_argument("-i2", "--rbp_path",
                        help= "path to RBP side data file (BED format)",
                        default="non")
    parser.add_argument("-r", "--list_of_replicats", action="store",
                        nargs='+',
                        dest="list_of_replicats", required=True,
                        help= "filenames list of all replicats")
    parser.add_argument("-o", "--out_path",
                        help= "path where output folder should be stored")
    parser.add_argument("-t", "--overlap_th",
                        help= "overlap threshold",
                        default="0.3")
    parser.add_argument("-s", "--score_th",
                        help= "score threshold",
                        default="0.5")
    parser.add_argument("-e", "--external_object",
                        help= "external rri  overlaping object (Interlap dict)",
                        default="non")
    parser.add_argument("-fh", "--filter_hybrind",
                        default="off",
                        help= "filter the data for hyprids alrady detected by ChiRA")
    parser.add_argument("-mo", "--mode",
                        default="train",
                        help= "function call within which cherri mode [train/eval]")



    args = parser.parse_args()
    input_path_RRIs = args.RRI_path
    file_rbp_pos = args.rbp_path
    replicats = args.list_of_replicats
    out_path = args.out_path
    overlap_th = args.overlap_th
    score_th = args.score_th
    external_object = args.external_object
    filter_hybrind = args.filter_hybrind
    mode = args.mode

    timestr = time.strftime("%Y%m%d")
    out_path =  out_path  + timestr + '_occ_out/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print('***added new folder***')

    # RRI thresholds
    # overlap_th = 0.3
    #score_th = 0.5
    # score_th = 'non'
    # RBP params
    seq_tag = '_RBP_side_'
    # context added to the T-> C side giving us the RBP interaction side
    context = 5
    exp_score_th = 10
    if file_rbp_pos == 'non':
        flag_prot = False
        print('no RBP occupyed postions given')
    else:
        flag_prot = True

    #### Get RRI data by calling find trusted RRI with a very low overlap th of 5%
    ### only take uniquly mapped reads but they do not need to be to stricke over the replicats:

    ####### Get RRI data
    rri_call_param = ('-i ' + input_path_RRIs + ' -r ' + ' '.join(replicats) +
                     ' -o ' + str(overlap_th) +' -n rri_occupied_regions -d ' +
                     out_path + ' -s ' +  str(score_th))
    if filter_hybrind == 'on':
        rri_call_param = f'{rri_call_param} -fh on'
    rri_call  = 'find_trusted_RRI.py '  + rri_call_param

    rri_file = (out_path + 'rri_occupied_regions_overlap_' +
                str(overlap_th) + '.csv')

    if len(replicats) == 1:
        print('Info: only one experiment is used to build occupyed regions')
        in_file = input_path_RRIs + replicats[0]
        print(in_file)
        # df_replicat = rl.read_chira_data(in_file)
        if mode == 'train':
            df_replicat = rl.read_chira_data(in_file, header='no', separater="\t")
        elif mode == 'eval':
            df_replicat = rl.read_chira_data(in_file, header='yes', separater=",")
        print(df_replicat)
        if score_th == 'non':
            df_rris = df_replicat
        else:
            df_filtered_replicat = rl.filter_score(df_replicat, float(score_th))
            #df_rris = rl.delet_empty_col(df_filtered_replicat)
            df_rris = df_filtered_replicat

            #print(df_rris.info())
    else:
        rl.call_script(rri_call)
        #print(rri_call)
        df_rris = rl.read_chira_data(rri_file, header='yes', separater=",")
        #out_path =  out_path_temp
    #df_rris = rl.read_chira_data(file_test, header='yes', separater=",")
    print(df_rris)
    inter_rep_two = build_interlap_occ_sides(df_rris, 'two')
    inter_rri = rl.mearge_overlaps(inter_rep_two, 'rri')

#check data:
    print('##RRI results ###')
    count_entrys(inter_rri, 'rri')

    ####### Get protein data
    if flag_prot:
        inter_rbp = get_prot_occ_regions(file_rbp_pos, exp_score_th, context)
        print('##RBP results ###')
        count_entrys(inter_rbp, 'rbp')

        # add the two inter laps together
        for key in inter_rri:
            if key in inter_rbp:
                inter_rri[key].add(list(inter_rbp[key]))

        #check data:
        print('##Results of both lists###')
        count_entrys(inter_rri, 'both')

    if external_object != 'non':
        inter_external = rl.load_occupyed_data(external_object)
        count_entrys(inter_external, 'external')
        for key in inter_rri:
            if key in inter_external:
                inter_rri[key].add(list(inter_external[key]))

    # save files
    or_path = out_path + "occupied_regions.obj"
    or_handle = open(or_path,"wb")
    pickle.dump(inter_rri,or_handle)
    or_handle.close()
    print('object contining InteLab object with start and end postions:\n%s'%or_path)

    # filter rri file and save:
    output_name = 'rri_occupied_regions' + '_overlapTH_' + str(overlap_th) + '_scoreTH_1.csv'
    # df_rris_filterd = df_rris[(df_rris.score_seq_1st_side >= 1) & (df_rris.score_seq_2end_side >= 1)]
    df_rris_filterd = rl.filter_score(df_rris, 1)

    df_rris_filterd = df_rris_filterd[df_rris_filterd['chrom_1st'] != False]
    df_final_output = df_rris_filterd[df_rris_filterd['chrom_2end'] != False]

    df_final_output.to_csv(out_path + output_name, index=False)



if __name__ == '__main__':
    main()
