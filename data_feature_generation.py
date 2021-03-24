#!/usr/bin/env python
import pandas as pd
import math
#import matplotlib as mpl
#import matplotlib.pyplot as plt
from collections import defaultdict
from interlap import InterLap
import sys
import argparse
import numpy as np
import rrieval.lib as rl


def read_data(in_file, header='no'):
    """
    Read RRI tabular file and convert to a dataframe including a header

        Parameters
        ----------
        in_file : tabular file with output of chira RRI results

        Raises
        ------
        nothing

        Returns
        -------
        df_interactions
            dataframe listing all interactions

        """
    df_temp = pd.read_table(in_file, header=None, sep="\t")
    # inclued header
    if header == 'no':
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
    # len(header)
    df_interactions = pd.DataFrame(df_temp.values, columns=header)
    return df_interactions

def filter_score(df_interactions):
    """
    Filter dataframe for instances with a score of 1

        Parameters
        ----------
        df_interactions : df including the containing all RRIs

        Raises
        ------
        nothing

        Returns
        -------
        df_interactions_single_mapped
            dataframes with instances filter for a score of 1

            """
    # filter input for score_seq_1st_side and score_seq_2end_side == 1
    df_interactions_single_mapped = df_interactions[(df_interactions.score_seq_1st_side == 1) & (df_interactions.score_seq_2end_side == 1)]
    #df_interactions_single_mapped
    return df_interactions_single_mapped


def get_chrom_list_no_numbers(df_interactions, chrom):
    """
    Generates a unique list of chromosmes or conticts for both interaction
    partners

        Parameters
        ----------
        df_interactions : df including the filtered RRIs
        chrom :  string indicating from wich seq the chromosome is

        Raises
        ------
        nothing

        Returns
        -------
        sort_list_chrom
            sorted list of unique chromosmes or contics which are not a number
            and present in the input data frame

        """
    chrom_list = df_interactions[chrom].unique().tolist()
    #convert all values to string in case it is not
    new_list = [str(el) for idx,el in enumerate(chrom_list)]
    sort_list_chrom = sorted(new_list)

    return sort_list_chrom


def get_list_chrom(df_interactions):
    """
    Generates a unique list of chromosmes or conticts for both interaction
    partners

        Parameters
        ----------
        df_interactions : df including the filtered RRIs

        Raises
        ------
        nothing

        Returns
        -------
        sort_list_chrom
            sorted list of unique chromosmes or contics which are not a number
            and present in the input data frame

        """
    chrom1_list = get_chrom_list_no_numbers(df_interactions, 'chrom_seq_1st_side')
    chrom2_list = get_chrom_list_no_numbers(df_interactions, 'chrom_seq_2end_side')
    list_chrom_no_int = list(set().union(chrom1_list,chrom2_list))
    sort_list_chrom = sorted(list_chrom_no_int)
    return sort_list_chrom


def build_interlap_for_replicat(df_interactions):
    """
    Building the inerlap objects for a fast overlap comparision for one replicat

        Parameters
        ----------
        df_interactions : df including the filtered RRIs

        Raises
        ------
        nothing

        Returns
        -------
        inter_rep
            inerlap objects for a fast overlap comparision

        """
    list_chrom_no_int = get_list_chrom(df_interactions)
    #print(list_chrom_no_int)
    first_key = ''
    # use defaultdict to key by chromosome.
    inter_rep = defaultdict(InterLap)

    for index, row in df_interactions.iterrows():
        row_content = row
        chrom_1 = str(row['chrom_seq_1st_side'])
        chrom_1_index = list_chrom_no_int.index(chrom_1)
        chrom_2 = str(row['chrom_seq_2end_side'])
        chrom_2_index = list_chrom_no_int.index(chrom_2)

        if chrom_2_index < chrom_1_index:
            first_key = str(chrom_2) + ':' + str(chrom_1)
            second_key = row['strand_seq_2end_side'] + ':' + row['strand_seq_1st_side']
            both_keys = first_key + ';' + second_key
            inter_rep[both_keys].add((row['start_seq_2end_side'], row['stop_seq_2end_side'], both_keys,
                                    [row['start_seq_1st_side'],row['stop_seq_1st_side']],
                                    ['swap',row['IntaRNA_prediction'], row['energy']], [row]))
        elif chrom_1_index <= chrom_2_index:
            first_key = str(chrom_1) + ':' + str(chrom_2)
            second_key = row['strand_seq_1st_side'] + ':' + row['strand_seq_2end_side']
            both_keys = first_key + ';' + second_key
            inter_rep[both_keys].add((row['start_seq_1st_side'], row['stop_seq_1st_side'], both_keys,
                                    [row['start_seq_2end_side'],row['stop_seq_2end_side']],
                                    ['no',row['IntaRNA_prediction'], row['energy']], row))
        else:
            print('error: something went wrong!!')

    return inter_rep




def build_replicat_library_to_compare(input_path, list_of_replicats):
    """
    Building for each replicat a inter object

        Parameters
        ----------
        input_path : path to the input files
        list_of_replicats: list of replicat file names

        Raises
        ------
        nothing

        Returns
        -------
        inter_replicat_list
            list for all replicat inter object
        no_replicats
            number of replicats

        """
    inter_replicat_list = []
    rep_size_list = []
    for file in list_of_replicats:
        in_file = input_path + '/' + file
        df_replicat = read_data(in_file)
        df_filtered_replicat = filter_score(df_replicat)
        rep_size = len(df_filtered_replicat)
        rep_size_list.append(rep_size)
        inter_replicat = build_interlap_for_replicat(df_filtered_replicat)
        inter_replicat_list.append(inter_replicat)
    #print(len(inter_replicat_list))
    # sort the replicat list by size...

    no_replicats = len(inter_replicat_list)
    inter_replicat_sorted_list = sort_list_replicat(inter_replicat_list, rep_size_list)

    return inter_replicat_sorted_list, no_replicats, rep_size_list


def sort_list_replicat(inter_replicat_list, rep_size_list):
    """
    looks for the rep with the least RRI and adds it to the first pair_position
    of the replicats list

        Parameters
        ----------
        inter_replicat_list : list contining all replicats

        Raises
        ------
        nothing

        Returns
        -------
        inter_replicat_sorted_list
            returns a list contaning the the inter opjects for each replicat

        """
    inter_replicat_sorted_list = inter_replicat_list
    min_rri_rep = min(rep_size_list)
    min_rri_rep_index = rep_size_list.index(min_rri_rep)
    inter_replicat_sorted_list[0], inter_replicat_sorted_list[min_rri_rep_index] = inter_replicat_list[min_rri_rep_index], inter_replicat_list[0]
    return inter_replicat_sorted_list


def rep_seq_pos(inter_rep):
    """
    extracts the start and stop positions of a interaction for one replicats
    for a inter opject:
    (s1, e1, 'Chrom1:Chrom2;strand1:strand2', [s2, e2], ['swap' or not, structure, enegy])

        Parameters
        ----------
        inter_rep :

        Raises
        ------
        nothing

        Returns
        -------
        rep_s1
            start postion of the first sequence
        rep_e1
            end postion of the first sequence
        rep_s2
            start postion of the second sequence
        rep_e2
            end postion of the second sequence

        """
    rep_s2 = inter_rep[3][0]
    rep_e2 = inter_rep[3][1]
    rep_s1 = inter_rep[0]
    rep_e1 = inter_rep[1]
    return rep_s1, rep_e1, rep_s2, rep_e2



def check_last_position(pos_replicat, inter_list, trusted_rri_temp,
                        pair_position, max_overlap):
    """
    Building for each replicat a inter object

        Parameters
        ----------
        pos_replicat: postion within the overlaping inter_inter list of the
            current overlapping repilcat
        inter_list: list of RRIs that overlap with with the first sequence
            of the first replicat
        trusted_rri_temp: list to store the rri of each replicat which has the
            highest overlap with the first replicat
        pair_position: postion of the replicat which has the highest overlap
            with the first replicat within the inter_list
        max_overlap: the overlap of the rri having the highest overlap among all
            overlapping rris between the first replicat and the current replicat


        Raises
        ------
        nothing

        Returns
        -------
        trusted_rri_temp
            list of rri that overlap including the one of the current replicat
            if on overlapping rri was found
        pair_position
            index postion of the rri is initilzed with -1 for the next
            replicat comparison round
        max_overlap
            overlap is initilzed with 0 here for the next replicat comparison
            round

        """
    if pos_replicat == len(inter_list):
        if pair_position >= 0:
            #print(inter_list)
            #print(pair_position)
            #print(inter_list[pair_position])
            trusted_rri_temp.append(inter_list[pair_position])
            #print('temp rri:')
            #print(trusted_rri_temp)
        pair_position = -1
        max_overlap = 0
    #print(trusted_rri_temp)
    return trusted_rri_temp, pair_position, max_overlap


def compare_overlap(overlap_pairs_list, overlap_th, no_replicats):
    """
    compare_overlap

        Parameters
        ----------
        overlap_pairs_list: list of overlapping RRIs
        overlap_th: pverlap threshold
        no_replicats: number of replicats

        Raises
        ------
        nothing

        Returns
        -------
        trusted_rri
            list contining the most overlaping rri for each replicat if there
            one could be found

        """
    max_overlap = 0
    pair_position = -1
    rep1 = overlap_pairs_list.pop(0)
    trusted_rri_temp = []
    trusted_rri_temp.append(rep1)
    #print(trusted_rri_temp)
    rep1_s1, rep1_e1, rep1_s2, rep1_e2 = rep_seq_pos(rep1)
    for inter_list in overlap_pairs_list:
        #print(inter_list)
        #print('beginn with %f and %f'% (max_overlap, pair_position))
        for idx,pair in enumerate(inter_list):
            #print(pair)
            #print('index:%i'%idx)
            s1, e1, s2, e2 = rep_seq_pos(pair)
            pos_replicat = idx+1
            #print(s1, e1, s2, e2)
            # Check if there is a overlap in the second sequence!
            if ((s2 <= rep1_s2 and e2 >= rep1_s2) or
                    (rep1_s2 <= s2 and rep1_e2 >= s2)):

                    # test if the overlap is big enaghe to continue!
                    overlap_rep1_rep2_seq1 = rl.calculate_overlap(rep1_s1,rep1_e1,s1,e1)
                    overlap_rep1_rep2_seq2 = rl.calculate_overlap(rep1_s2,rep1_e2,s2,e2)
                    #print(overlap_rep1_rep2_seq2, overlap_rep1_rep2_seq1)
                    #print(overlap_th)
                    if ((overlap_rep1_rep2_seq2 < overlap_th) or (overlap_rep1_rep2_seq1 < overlap_th)):
                        # did not full fill the needed overlap th
                        #print('not enagh overlap!')
                        trusted_rri_temp, pair_position, max_overlap = check_last_position(pos_replicat,
                                                                                           inter_list,
                                                                                           trusted_rri_temp,
                                                                                           pair_position,
                                                                                           max_overlap)
                    else:
                        mean_overlap = (overlap_rep1_rep2_seq1+overlap_rep1_rep2_seq2)/2
                        #print(mean_overlap)
                        if max_overlap < mean_overlap:
                            max_overlap = mean_overlap
                            #print('max overlap: %f'%max_overlap)
                            pair_position = idx
                            #print('at positon %i'%idx)
                        #print('length of list: %i|| index: %i'%(len(inter_list), (idx+1))
                        #print('length of list: %i|| index: %i'%(len(inter_list), pos_replicat)
                        #print(pos_replicat,len(inter_list),pair_position)
                        #print('enagh overlap!!')
                        trusted_rri_temp, pair_position, max_overlap = check_last_position(pos_replicat,
                                                                                           inter_list,
                                                                                           trusted_rri_temp,
                                                                                           pair_position,
                                                                                           max_overlap)
            else:
                # in case the last entry of the list is not overlapping!
                #print('not overlapping at all!')
                #print('start and stop or rep %i and %i'% (rep1_s2, rep1_e2))
                trusted_rri_temp, pair_position, max_overlap = check_last_position(pos_replicat,
                                                                                   inter_list,
                                                                                   trusted_rri_temp,
                                                                                   pair_position,
                                                                                   max_overlap)
    if len(trusted_rri_temp) == no_replicats:
        trusted_rri = trusted_rri_temp
    else:
        trusted_rri = ['nan']
    #print(trusted_rri_temp)
    #print(trusted_rri)
    return trusted_rri



def find_relayble_replicats(inter_replicat_list, overlap_th, no_replicats):
    """
    find_relayble_replicats

        Parameters
        ----------
        inter_replicat_list: list of replicat file names
        overlap_th : overlap threshold
        no_replicats: number of replicats

        Raises
        ------
        nothing

        Returns
        -------
        no_relayble_rri
            number of relayble RRIs found in all replicats
        trusted_rri_list
            list of relayble RRIs found in all replicats
        no_replicats
            number of replicats

        """
    # get the first replicat
    # no_replicats = len(inter_replicat_list)
    inter_rep1 = inter_replicat_list.pop(0)
    overlap_pairs_list = []
    trusted_rri_list = []
    no_relayble_rri = 0
    #print(no_replicats)

    # find the pairs in all replicats
    for key in inter_rep1:
        #print(inter_rep1[key])
        for rep1 in inter_rep1[key]:
            #print(rep1)
            overlap_pairs_list = []
            overlap_pairs_list.append(rep1)
            # now compare to other rep
            counter = 1
            for idx, inter_rep_next in enumerate(inter_replicat_list):
                #print(list(inter_rep_next[key].find(rep1)))
                #inter_rep_next[key]
                if (rep1 in inter_rep_next[key]):
                    # append the RRI pairs found in the next replicat
                    #print(inter_rep_next[key])
                    overlap_pairs_list.append(list(inter_rep_next[key].find(rep1)))
                    #print(idx)
                    if (idx == (no_replicats -2)):
                        overlap_pairs_list
                        # check if all overlaps in overlap_pairs_list are in overlap_th
                        trusted_rri = compare_overlap(overlap_pairs_list, overlap_th, no_replicats)
                        #print('final rri list:')
                        #print(trusted_rri)
                        if len(trusted_rri) <= 1:
                            continue
                        elif len(trusted_rri) == no_replicats:
                            no_relayble_rri += 1
                            trusted_rri_list.append(trusted_rri)
                            #print('trusted RRI:')
                            #print(trusted_rri_list)
                        else:
                            print('error something went wrong!!!')



                else:
                    # we do not need to evaluate this RRI instance further, because it does not overlap
                    break
    return no_relayble_rri, trusted_rri_list, no_replicats


def get_numbers_nan(no_replicats, trusted_rri_list):
    """
    finding the interacitons, where the enegry value for some or all
    interaction could not be computed. Mening no IntaRNA prediciton is
    avalable.

        Parameters
        ----------
        overlap_th : overlap threshold
        list_of_replicats: list of replicat file names

        Raises
        ------
        nothing

        Returns
        -------
        instances_just_nan_list
            list of all interactions without a IntaRNA perdiction
        instances_also_nan_list
            list of all interactions missing the IntaRNA perdiction for some
            replicats
        instances_no_nan_list
            list of all interactions where all replicats have a
            IntaRNA perdiction

        """
    instances_just_nan_list = []
    instances_also_nan_list = []
    instances_no_nan_list = []
    for rep_list in trusted_rri_list:
        sturcture_list = []
        count_nan = 0
        for rep_instans in rep_list:
            structure = rep_instans[4][1]
            #print(type(structure))
            if str(structure) == 'nan':
                count_nan +=1
            sturcture_list.append(structure)
        if count_nan == 0:
            instances_no_nan_list.append(rep_list)
        elif count_nan == no_replicats:
            instances_just_nan_list.append(rep_list)
        else:
            instances_also_nan_list.append(rep_list)

    return instances_just_nan_list, instances_also_nan_list, instances_no_nan_list


def get_seq_IntaRNA_calls(trusted_rri_list):
    """
    collect all

        Parameters
        ----------
        trusted_rri_list : list of trusted rris! each rri is a list of
        [inter_instance rep1, [list overlapping second rep]...,
        [list overlapping last rep]]

        Raises
        ------
        nothing

        Returns
        -------
        instances_just_nan_list
        """

    nan_seq_list = []
    for rep_list in trusted_rri_list:
        for rep_instans in rep_list:
            enegy = rep_instans[4][2]
            #print(enegy)
            if str(enegy) == 'nan':
                #print(rep_instans[5][0])
                if isinstance(rep_instans[5][0], int):
                    #print('\nInterger: %i\n'%rep_instans[5][0])
                    nan_seq_list.append(rep_instans[5])
                else:
                    nan_seq_list.append(rep_instans[5][0])
                    #print('\nshould be series:')
                    #print(rep_instans[5][0])
    df_output = concat_series_objects(nan_seq_list)

    return df_output


def get_enegy_seqlen(trusted_rri_list):
    """
    finding the interacitons, where the enegry value for some or all
    interaction could not be computed. Mening no IntaRNA prediciton is
    avalable.

        Parameters
        ----------
        trusted_rri_list : list of trusted rris! each rri is a list of
        [inter_instance rep1, [list overlapping second rep]...,
        [list overlapping last rep]]

        Raises
        ------
        nothing

        Returns
        -------
        instances_just_nan_list
        """

    enegy_list = []
    interaction_length = []
    final_output_list = []

    for rep_list in trusted_rri_list:
        enegy_temp_list = []

        #print(rep_list)
        for rep_instans in rep_list:
            enegy = rep_instans[4][2]
            if str(enegy) == 'nan':
                enegy_temp_list.append(100)
            else:
                enegy_temp_list.append(enegy)
        min_enegry = min(enegy_temp_list)
        #print(enegy_temp_list)
        min_enegry_index = enegy_temp_list.index(min_enegry)
        min_enegy_rep = rep_list[min_enegry_index]
        #print(min_enegy_rep)
        #print(min_enegy_rep[5])
        instance = min_enegy_rep[5][0]
        #print(instance)
        if isinstance(min_enegy_rep[5][0], int):
            #print('Instance is int: %i'%instance)
            #rint(min_enegy_rep[5])
            instance = min_enegy_rep[5]
        else:
            instance = min_enegy_rep[5][0]
        #print('Instance element\n: %s\n'%instance)
        final_output_list.append(instance)
        interaction_length = get_seq_lengths(min_enegy_rep, interaction_length)
        enegy_list.append(min_enegy_rep[4][2])

    df_output = concat_series_objects(final_output_list)
    #print(df_output.info())
    return enegy_list, interaction_length, df_output


def get_seq_lengths(min_enegy_rep, interaction_length):
    """
    concatinate all series instances with on list. Since the series objects are
    the content of one row the dataframe is transposed
        Parameters
        ----------
        min_enegy_rep : Inter instance containing all sequence positons
        interaction_length: list with length of the two interaction RNAs so far

        Raises
        ------
        nothing

        Returns
        -------
        interaction_length:
            list where the length of the two interaction RNAs are appendend

        """
    s1, e1, s2, e2 = rep_seq_pos(min_enegy_rep)
    sequence1_len = sequence_length(s1, e1)
    sequence2_len = sequence_length(s2, e2)
    interaction_length.append((sequence1_len, sequence2_len))
    return interaction_length




def concat_series_objects(list_of_series):
    """
    concatinate all series instances with on list. Since the series objects are
    the content of one row the dataframe is transposed
        Parameters
        ----------
        list_of_series : list of seris objects containing the row information
        of the choosen trusted RRI

        Raises
        ------
        nothing

        Returns
        -------
        df_output:
            dataframe containg all

    >>> s1 = pd.Series([1, 2], index=['A', 'B'], name='s1')
    >>> s2 = pd.Series([3, 4], index=['A', 'B'], name='s2')
    >>> concat_series_objects([s1, s2])

        """
    #for i in list_of_series:
        #print(i[0])
    df_temp = pd.concat(list_of_series, axis=1)
    df_output = df_temp.T
    return df_output




def sequence_length(start, end):
    """
    compute sequence length

        Parameters
        ----------
        start: start position of a sequence
        end: end postion of a sequence

        Raises
        ------
        nothing

        Returns
        -------
        seq_length
            length of a sequence


        """
    seq_length = end - start +1

    return seq_length





def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_path", action="store", dest="input_path",
                        required=True,
                        help= "path to folder storing all input data")
    parser.add_argument("-r", "--list_of_replicats", action="store",
                        nargs='+',
                        dest="list_of_replicats", required=True,
                        help= "list having filenames of all replicats")
    parser.add_argument("-o", "--overlap_th", action="store",  type=float,
                        dest="overlap_th", required=True,
                        help= "overlap threshold to find relyable RRIs")
    parser.add_argument("-d", "--output_path", action="store", dest="output_path",
                        required=True,
                        help= "path output reposetory")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of positve trusted RRIs")


    args = parser.parse_args()
    input_path = args.input_path
    list_of_replicats = args.list_of_replicats
    #print(list_of_replicats)
    overlap_th = args.overlap_th
    output_path = args.output_path
    experiment_name = args.experiment_name


    plot_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/RNA_RNA_binding_evaluation/plots/'
    #output_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/'
    #output_tag = 'PARIS_mES'
    output_name = experiment_name + '_overlap_' + str(overlap_th) + '.cvs'

    plot_path_full = plot_path + '_' + str(overlap_th) + '_'


    inter_replicat_list, no_replicats, rep_size_list = build_replicat_library_to_compare(input_path, list_of_replicats)
    len_smalles_replicat = rep_size_list[0]
    print(rep_size_list)
    print(len_smalles_replicat)
    no_relayble_rri, trusted_rri_list, no_replicats = find_relayble_replicats(inter_replicat_list, overlap_th, no_replicats)
    instances_just_nan_list, instances_also_nan_list, instances_no_nan_list = get_numbers_nan(no_replicats, trusted_rri_list)
    enegy_list, interaction_length, df_output_temp = get_enegy_seqlen(instances_no_nan_list)
    enegy_list_also_nan, interaction_length_also_nan, df_output_nan_temp = get_enegy_seqlen(instances_also_nan_list)
    out_for_IntaRNA_calls_df = get_seq_IntaRNA_calls(instances_also_nan_list)
    #print(out_for_IntaRNA_calls)

    # print(df_output_temp.info())
    # print(df_output_nan_temp.info())

    df_final_output = pd.concat([df_output_temp, df_output_nan_temp])
    df_final_output.to_csv(output_path + output_name, index=False)
    out_for_IntaRNA_calls_df.to_csv(output_path + 'NAN_' + output_name, index=False)

    percentage_trustable_rri_all = no_relayble_rri/len_smalles_replicat


    print('######\n for %i replicates the following number of reliable interactions are found: %i (%f)'%(no_replicats, no_relayble_rri, percentage_trustable_rri_all))
    print('the distribution of the interactions are:')
    print('Number of only nan RRI: %i'%len(instances_just_nan_list))
    print('Number of some nan RRI: %i'%len(instances_also_nan_list))
    print('Number of no nan RRI: %i'%len(instances_no_nan_list))
    #print(no_relayble_rri)
    #print(len_smalles_replicat)
    print('######')



    #### Plotting ######

    #histogrom enegy
    #fig1 = plt.figure()
    #bins = np.arange(min(enegy_list), max(enegy_list), 5)

    #plt.hist(enegy_list, bins=bins)
    #fig1.savefig(plot_path + "histogram_enegy.pdf", bbox_inches='tight')

    #seq1_len_list = [len[0] for len in interaction_length_also_nan]
    #seq2_len_list = [len[1] for len in interaction_length_also_nan]

    #d = {'rri_seq1': seq1_len_list, 'rri_seq2': seq2_len_list}
    #df_rri_len = pd.DataFrame(data=d)

    #myFig = plt.figure()
    #boxplot = df_rri_len.boxplot(column=['rri_seq1', 'rri_seq2'])

    #myFig.savefig(plot_path + "boxplot_rri_len_seq.pdf", bbox_inches='tight')

    # input_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/RNA_RNA_binding_evaluation/data/training/Paris/'
    # list_of_replicats = ['test_rep1.tabular', 'test_rep2.tabular', 'test_rep3.tabular']

if __name__ == '__main__':
    main()