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


        Returns
        -------
        inter_rep
            inerlap objects for a fast overlap comparision

        """

    first_key = ''
    flag_name = 'crl'
    # use defaultdict to key by chromosome.
    inter_rep = defaultdict(InterLap)

    if flag_name == 'hybrid':
        names_first = ['chrom_seq_1st_side', 'strand_seq_1st_side', 'start_seq_1st_side', 'stop_seq_1st_side']
        names_second = ['chrom_seq_2end_side', 'strand_seq_2end_side', 'start_seq_2end_side', 'stop_seq_2end_side']
    elif flag_name == 'crl':
        names_first = ['chrom_1st','strand_1st','start_1st','end_1st' ]
        names_second = ['chrom_2end','strand_2end','start_2end','end_2end']

    if flag == 'two':
        list_chrom_no_int = rl.get_list_chrom(df_interactions)
        df_interactions[names_first[0]] = df_interactions[names_first[0]].apply(lambda x: rl.check_convert_chr_id(x))
        df_interactions[names_second[0]] = df_interactions[names_second[0]].apply(lambda x: rl.check_convert_chr_id(x))
        for index, row in df_interactions.iterrows():
            row_content = row
            if not row[names_first[0]]:
                print('can not use chromosome')
            else:
                both_keys1 = str(row[names_first[0]]) + ';' + row[names_first[1]]
                inter_rep[both_keys1].add((row[names_first[2]], row[names_first[3]], [row]))

            if not row[names_second[0]]:
                print('can not use chromosome')
            else:
                both_keys2 = str(row[names_second[0]]) + ';' + row[names_second[1]]
                inter_rep[both_keys2].add((row[names_second[2]], row[names_second[3]], [row]))


    elif flag == 'one':
        list_chrom_no_int = rl.get_chrom_list_no_numbers(df_interactions, 'chrom')
        for index, row in df_interactions.iterrows():
            row_content = row

            both_keys = str(row['chrom']) + ';' + row['strand']
            inter_rep[both_keys].add((row['start'], row['end'], [row]))

    return inter_rep





def join_pos(pos_list):
    """
    join positons will join start end end postions whick are overlaping

        Parameters
        ----------
        pos_list : list of tupels containg (start, end) position
        info: information what inter object

        Returns
        -------
        inter_obj_new
            inerlap objects with the mearged positons

    >>> join_pos([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]
    >>> join_pos([(2, 6), (4, 10)])
    [(2, 10)]
    !!! COPY !!!
    """
    if len(pos_list) < 2: return pos_list
    pos_list.sort()
    joint_pos_list = [pos_list[0]]
    for next_i, (s, e) in enumerate(pos_list, start=1):
        if next_i == len(pos_list):
            joint_pos_list[-1] = joint_pos_list[-1][0], max(joint_pos_list[-1][1], e)
            break

        ns, ne = pos_list[next_i]
        if e > ns or joint_pos_list[-1][1] > ns:
            joint_pos_list[-1] = joint_pos_list[-1][0], max(e, ne, joint_pos_list[-1][1])
        else:
            joint_pos_list.append((ns, ne))
    return joint_pos_list


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--RRI_path",
                        help= "path to folder storing all RRI data (tabel)",
                        default="/vol/scratch/data/RRIs/Paris/")
    parser.add_argument("-i2", "--rbp_path",
                        help= "path to RBP side data file (bed format)",
                        default="/vol/scratch/data/human_RBP_coverage/GSE38355_ProtOccProf_4SU_consensus_TC_hg38.bed")
    parser.add_argument("-r", "--list_of_replicats", action="store",
                        nargs='+',
                        dest="list_of_replicats", required=True,
                        help= "list having filenames of all replicats")
    parser.add_argument("-o", "--out_path",
                        help= "path to folder storing outputfiles",
                        default="/vol/scratch/data/RRIs/")



    args = parser.parse_args()
    input_path_RRIs = args.RRI_path
    file_rbp_pos = args.rbp_path
    replicats = args.list_of_replicats
    out_path = args. out_path

    timestr = time.strftime("%Y%m%d")
    out_path = '/' + out_path + '/' + timestr + '_occ_out/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # RRI thresholds
    overlap_th = 0.3
    score_th = 0.5
    # RBP params
    seq_tag = '_RBP_side_'
    # context added to the T-> C side giving us the RBP interaction side
    context = 20
    exp_score_th = 10

    #### Get RRI data by calling find trusted RRI with a very low overlap th of 5%
    ### only take uniquly mapped reads but they do not need to be to stricke over the replicats:

    ####### Get RRI data
    rri_call_param = ('-i ' + input_path_RRIs + ' -r ' + ' '.join(replicats) + ' -o ' +
                     str(overlap_th) +' -n rri_occupied_regions -d ' + out_path+
                     ' -s ' +  str(score_th))
    rri_call  = 'python -W ignore find_trusted_RRI.py '  + rri_call_param

    rri_file = (out_path + 'rri_occupied_regions_overlap_' +
                str(overlap_th) + '.cvs')

    rl.call_script(rri_call)

    #print(rri_call)
    #file_test =  '/vol/scratch/data/trusted_RRIs/paris_HEK293T_overlap_0.6_snoRNA.cvs'
    df_rris = rl.read_chira_data(rri_file, header='yes', separater=",")
    #df_rris = rl.read_chira_data(file_test, header='yes', separater=",")
    #print(df_rris)
    inter_rep_two = build_interlap_occ_sides(df_rris, 'two')


    #for key in inter_rep_two:
        #print(key)
        #print(list(inter_rep_two[key]))

    ####### Get protein data

    ## Read protein data in bed format
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
    #print(df_context)

    inter_rep_one = build_interlap_occ_sides(df_context, 'one')
    # print('$$$$$$$$$$')
    # for key in inter_rep_one:
    #     print(key)
    #     print([(i[0],i[1]) for i in list(inter_rep_one[key])])
    # print('$$$$$$$$$$')
    inter_rri = rl.mearge_overlaps(inter_rep_two, 'rri')
    inter_rbp = rl.mearge_overlaps(inter_rep_one, 'rbp')

#check data:
    print('##RRI results ###')
    count2 = 0
    for key in inter_rep_two:
        count2 += len(list(inter_rep_two[key]))
        #print(key)
        #print([(i[0],i[1]) for i in list(inter_rep_two[key])])
    print('##########')
    print('entrys in List rri: %i' %count2)

    count3 = 0
    for key in inter_rbp:
        count3 += len(list(inter_rbp[key]))
        # print(key)
        # print(list(inter_rbp[key]))
    print('##########')
    print('entrys in List rbp: %i' %count3)



    # add the two inter laps together
    for key in inter_rri:
        if key in inter_rbp:
            inter_rri[key].add(list(inter_rbp[key]))

    #check data:
    count = 0
    for key in inter_rri:
        count += len(list(inter_rri[key]))
        # print(key)
        # print(list(inter_rri[key]))
    print('##########')
    print('entrys in List both: %i' %count)

    # save files
    or_path = out_path + "/occupied_regions.obj"
    or_handle = open(or_path,"wb")
    pickle.dump(inter_rri,or_handle)
    or_handle.close()
    print('object contining InteLab object with start and end postions:\n%s'%or_path)



if __name__ == '__main__':
    main()
