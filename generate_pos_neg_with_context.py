#!/usr/bin/env python
import pandas as pd
import math
import numpy as np
#import matplotlib as mpl
from collections import defaultdict
from interlap import InterLap
import sys
from operator import itemgetter
#import seaborn
import argparse
#import csv
#import collections
#from operator import itemgetter
import subprocess
import random
import re
import os
import time
import rrieval.lib as rl
import pickle

def check_context_extention(df, context, output_path, context_not_full):
    """
    assert if the full context was added for all sequences.

        Parameters
        ----------
        df: dataframe holding the extended context sequences
        context: amount of nt added on both sides
        output_path:

        Raises
        ------
        nothing

        """
    # check if context is fully added:
    df['seq_len'] = df['ineraction_side_1st'].astype(str).map(len)
    df['seq_con_len'] = df['con_target'].astype(str).map(len)
    #df_test = df[(df['seq_len']+(2*context)) != (df['seq_con_len'].all())]
    df['con_target'] = df['con_target'].replace('', np.nan)
    df['con_query'] = df['con_query'].replace('', np.nan)
    #null_col_target = df['con_seq_only_target'].isnull().sum()
    #null_col_query = df['con_seq_only_query'].isnull().sum()

    df_filtered_target = df.dropna(axis=0, subset=['con_target'])
    df_filtered_query = df_filtered_target.dropna(axis=0, subset=['con_query'])
    # print(len(df))
    # print(len(df_filtered_query))
    null_col = len(df)-len(df_filtered_query)

    if null_col > 0:
        print('Warning: context was not fully added for %i RRIs' %(null_col))
        context_not_full += 1
    return df_filtered_query, context_not_full


def get_context_pos(df, context, start, end, name):
    """
    get find object for Interlab by adding the context to start and end postions

        Parameters
        ----------
        df: dataframe holding the extended context sequences
        context: amount of nt added on both sides

        Returen
        ------
        df with a added coulem find_target and find_query

        """
    col_name_s = name + '_con_s'
    col_name_e = name + '_con_e'

    df[col_name_s] = df[start] - context
    df[col_name_e] = df[end] + context

    return df


def extention_df(df):
    """
    defining colum with ID and empty colums to store the context sequences

        Parameters
        ----------
        df: dataframe

        Raises
        ------
        nothing

        Returns
        -------
        df
            colum update dataframe

        """
    # add RRI number as ID
    df['interaction_no'] = np.arange(len(df))
    #add ids to the df
    df['ID1']=  df['chrom_1st'].astype(str) + ':' + df['start_1st'].astype(str)+ ':' + df['end_1st'].astype(str)+ ':' + df['strand_1st'].astype(str)+ ':' + df['interaction_no'].astype(str)
    df['ID2']=  df['chrom_2end'].astype(str) + ':' + df['start_2end'].astype(str)+ ':' + df['end_2end'].astype(str)+ ':' + df['strand_2end'].astype(str)+ ':' + df['interaction_no'].astype(str)
    # add coulums for the sequeces inculding the context
    df['con_target'] = ''
    df['con_query'] = ''
    return df


def find_occu_overlaps(dict_key,s,e, occupyed_InteLab):
    """
    finding all interaction sides for a given sequence by genomic indexes

        Parameters
        ----------
        dict_key: key chrom:strand of the given sequence neede for the occupyed_InteLab
        s: start postion of given sequence
        e: end postion of given sequence
        occupyed_InteLab: Interlap object having all interacting side postions

        Raises
        ------
        nothing

        Returns
        -------
        occupied_regions_list
            list contining all occupied positions for the given sequence

        """

    temp_list = list(occupyed_InteLab[dict_key].find((s,e)))
    #print('###Overlaping list ########')
    #print(dict_key)
    #print(s,e)
    #print(list(occupyed_InteLab[dict_key]))
    if len(temp_list) < 1:
        print('Warning: occupyed list is not compleate, please repuduce the occupied object again.')
        sys.exit()
    if temp_list:
        #print('found overlap:')
        occupied_regions_list = temp_list
        # occupied_regions_list.append(temp_list)
        #print(occupied_regions_list)
    return occupied_regions_list


def convert_positions(s_seq, e_seq, s_side, e_side, strand):
    """
    return postions one based of side in relation to given sequence

        Parameters
        ----------
        s_seq: start postion of given sequence
        e_seq: end postion of given sequence
        s_side: start postion of given occupied side
        e_side: start postion of given occupied side
        strand: + or -

        Returns
        -------
        new_e_side
            one based end postion of occupied side
        new_s_side
            one based start postion of occupied side
     >>> convert_positions(20,30,22,24,'+')
     [3, 4]
     >>> convert_positions(20,30,12,34,'+')
     Warning: full context is overlaped
     [1, 10]
     >>> convert_positions(20,30,22,24,'+')
     [3, 4]
        """
    # check than side are outside the sequence
    if (s_side <= s_seq) or (e_side >= e_seq):
        if (s_side <= s_seq) and (e_side >= e_seq):
            print('Warning: full context is overlaped')
            new_end = e_seq - s_seq
            return  [1, new_end]
        # changing start
        elif s_side <= s_seq:
            s_side = s_seq
            if e_side == s_side:
                print('Warning: start and end is same positon START!!!')
                # should break!!
                return(['nan','nan'])
        elif e_side > e_seq:
            e_side = e_seq
    # compute postions
    if strand == '+':
        # new_end = e_seq - s_seq
        if e_side == s_side:
            print('Warning: start and end is same positon!!!+END')
            return[(e_side-s_seq),(e_side-s_seq)]
        new_e_side = e_side - s_seq
        new_s_side = (s_side+1) - s_seq
    elif strand == '-':
        if e_side == s_side:
            print('Warning: start and end is same positon!!!-END')
            return[1,1]
        new_e_side = e_seq - (s_side+1) +1
        new_s_side = e_seq - e_side +1
    else:
        print('strand was not found correctly')
    return [new_s_side, new_e_side]


def decode_Intarna_output(out):
    """
    Intarna_call

        Parameters
        ----------
        out: IntaRNA terminal output


        Returns
        -------
        df
            df incuding IntaRNA result

        """

    # out, err = process.communicate()
    out = out.decode('utf-8').strip().split('\n')
    #print(out)
    for idx, line in enumerate(out):
        #print(idx)
        line = line.strip().split(';')
        #line = line.strip()
        #line = line.split(';')
        #print(line)
        if idx == 0:
            col = line
            result = 'empty'
        elif idx == 1:
            result = 'exists'
            df = pd.DataFrame([line], columns=col)
        elif idx > 1:
            df_two = pd.DataFrame([line], columns=col)
            df = pd.concat([df, df_two])
        #print(line)
    #print(result)
    if result == 'empty':
        no_values = ['nan']*len(col)
        df = pd.DataFrame([no_values], columns=col)

    return df


def get_neg_pos_intarna_str(occupied_regions, neg_param, pos_s, pos_e):
    """
    get_neg_pos_intarna_st

        Parameters
        ----------
        occupied_regions:
        neg_param:
        pos_s:
        pos_e:


        Returns
        -------
        neg_param
            string for IntaRNA call

        """
    for idx, i in enumerate(occupied_regions):

        s = i[0]
        e = i[1]
        if (s == pos_s and e == pos_e and idx == 0):
            neg_param = ''
            break
        elif idx == 0:
            positions = str(i[0]) + '-' + str(i[1])
            neg_param = neg_param + positions
        elif (s == pos_s and e == pos_e):
            continue
        else:
            positions = str(i[0]) + '-' + str(i[1])
            neg_param = neg_param + ',' + positions
    neg_param = neg_param + '\"'
    return neg_param



def join_result_and_infos(df, lost_inst, row, list_rows_add):
    """
    The function checks whether IntaRNA could predict anything. If it did
    than additional infomation are appende to the inarna results staming from
    the original trusted RRI instace. This additional coulum naes are stored
    in the list_rows_add list.

        Parameters
        ----------
        df: IntaRNA result df
        lost_inst: no of already lost instances
        row: row of RRI dataframe
        list_rows_add: header line for to extract from row


        Returns
        -------
        df_result
            results of row and IntaRNA def
        lost_inst
            number of instances lost
        """

    if df['hybridDP'][0] == 'nan':
        lost_inst += 1
        return df, lost_inst
    else:
        val = {}
        # appeding information for the RRI that is not predicted by IntaRNA
        for col_n in list_rows_add:
            #print(col_n)
            val[col_n] = [row[col_n]]*len(df)

        df_append = pd.DataFrame.from_dict(val)
    #print(df_append)
    df_result = pd.merge(df, df_append, left_index=True, right_index=True)
    return df_result, lost_inst


def add_block_end_pos(seq_s, seq_e, block_ends,strand, new_occ_list):
    """
    convert the list conining all occupied postions

        Parameters
        ----------
        seq_s: target/query sequence start pos
        seq_e: target/query sequence end pos
        block_ends: nt that should be blocked at the end and start of seq
        strand: +/-
        new_occ_list: list with occupyed sequences

        Returns
        -------
        new_occ_list
            list start and stop postions for the end postions
        """
    block_ends_s = seq_s + block_ends
    block_ends_e = seq_e - block_ends
    block_front = convert_positions(seq_s, seq_e, seq_s, block_ends_s, strand)
    new_occ_list.append(block_front)
    block_back = convert_positions(seq_s, seq_e, block_ends_e, seq_e, strand)
    new_occ_list.append(block_back)

    return new_occ_list



def get_pos_occ_list(seq_s, seq_e, seed_s, seed_e, strand, converted_occ_list):
    """
    convert the list conining all occupied postions

        Parameters
        ----------
        seq_s: target/query sequence start pos
        seq_e: target/query sequence end pos
        seed_s: target/query seed sequence start pos
        seed_e: target/query seed sequence end pos
        strand: +/-
        converted_occ_list: list with occupyed sequences

        Returns
        -------
        pos_occ_list
            list not containing seed postions
        >>> get_pos_occ_list(1, 100, 50, 59, '+', [(2,10),(55,60),(62,70)])
        [(2, 10), (62, 70)]
        >>> get_pos_occ_list(1, 100, 30, 40, '+', [(30,40)])
        []
        >>> get_pos_occ_list(1, 100, 30, 40, '+', [(25,30), (39,45), (60,70)])
        [(60, 70)]
        """
    seed_pos = convert_positions(seq_s, seq_e, seed_s, seed_e, strand)
    pos_occ_list = []
    #print(converted_occ_list)
    for positions in converted_occ_list:
        s = positions[0]
        e = positions[1]
        if s >= seed_pos[0] and s <= seed_pos[1]:
            # found seed overlap
            continue
        elif e >= seed_pos[0] and e <= seed_pos[1]:
            # found seed overlap
            continue
        else:
            pos_occ_list.append(positions)

    return pos_occ_list




def convert_occu_positons(seq_s, seq_e, occupied_regions,strand, block_ends, seed_s, seed_e):
    """
    convert the list conining all occupied postions

        Parameters
        ----------
        seq_s: target/query sequence start pos
        seq_e: target/query sequence end pos
        occupied_regions: occupied regions list [(s,e)(s,e)...]
        strand: +/-
        block_ends: nt that should be blocked at the end and start of seq

        Returns
        -------
        new_occ_list
            list with updated postions
        """
    new_occ_list = []
    for i in occupied_regions:
        # print(i)
        s = i[0]
        e = i[1]
        pos_new = convert_positions(seq_s, seq_e, s, e, strand)
        if pos_new[0] == 'nan':
            continue
        else:
            new_occ_list.append(pos_new)

    if block_ends > 0:
        new_occ_list = add_block_end_pos(seq_s, seq_e, block_ends,strand,
                                             new_occ_list)

    new_occ_list = sorted(new_occ_list, key=itemgetter(0))
    occ_list_temp = [(i[0],i[1]) for i in new_occ_list]
    #print(inter_list_temp)
    mearged_list = rl.join_pos(occ_list_temp)
    #print(inter)
    converted_occ_list = [(i[0],i[1], 'converted') for i in mearged_list]

    pos_occ_list = get_pos_occ_list(seq_s, seq_e, seed_s, seed_e,
                                    strand, converted_occ_list)

    return converted_occ_list, pos_occ_list



def decode_IntaRNA_call(call, lost_inst, row, list_rows_add, df_data, no_sub_opt, no_less_sub_opt):
    """
    decode the IntaRNA call

        Parameters
        ----------
        call: IntaRNA call
        lost_inst: no of lost instances so fare
        row: infos of the curren RRI (line of df)
        list_rows_add:
        df_data: df contining the already calculated training instances
        no_sub_opt: suboptimals that should reported by IntaRNA if possible
        no_less_sub_opt: no of suboptimals that could not be predicted

        Returns
        -------
        df_data
            calculated instaces including the new prediction if prediction was sucessfull
        lost_inst
            number of lost instaces so fare. Inceased of one if prediction failed
        no_less_sub_opt
            no of suboptimal that could not be predicted so fare!
        """
    out = rl.call_script(call,reprot_stdout=True)
    #print(call)
    df = decode_Intarna_output(out)
    #print(df)
    # reset the indix of the df and avoid the old index being added as a column
    df = df.reset_index(drop=True)

    df_result, lost_inst_new = join_result_and_infos(df,
                                                         lost_inst,
                                                         row, list_rows_add)
    if lost_inst < lost_inst_new:
        # IntaRNA could not predict anything!
        lost_inst = lost_inst_new
    else:
        #print('instace appended to data')
        # check if found all number of subotpimals
        if no_sub_opt != len(df_result):
            print('Warning IntaRNA could not find %i interaction but %i interactions for a particular call'%(no_sub_opt, len(df_result)))
            no_less_sub_opt += no_sub_opt - len(df_result)

        df_data = pd.concat([df_data, df_result])
    return df_data, lost_inst, no_less_sub_opt


def get_context_added(input_rris, output_path, genome_file, context, context_not_full,context_file,chrom_len_file):
    """
    get_context

        Parameters
        ----------
        input_rris:
        output_path:
        genome_file:
        context:
        context_not_full:
        context_file:

        Returns
        -------
        df_contex
            datafram including the context appended sequences
        """
    df_RRIs = pd.read_table(input_rris, sep=",")
    print(len(df_RRIs))

    # adding context by including infors into the df
    df_RRIs = extention_df(df_RRIs)
    df_target = rl.get_context('target', df_RRIs, output_path,
                                genome_file, context, chrom_len_file)
        #print(df_target)
    df_context_seq = rl.get_context('query', df_target, output_path,
                                    genome_file, context, chrom_len_file)
        # print(df_context)


    df_filted_RRIs, context_not_full = check_context_extention(df_context_seq,
                                                               context,
                                                               output_path,
                                                               context_not_full)
    #### context added df saved!

    df_contex = get_context_pos(df_filted_RRIs, context, 'start_1st',
                                    'end_1st', 'target')
    df_contex = get_context_pos(df_contex, context, 'start_2end',
                                    'end_2end', 'query')

    df_contex['chrom_1st'] = df_contex['chrom_1st'].apply(lambda x: rl.check_convert_chr_id(x))
    df_contex['chrom_2end'] = df_contex['chrom_2end'].apply(lambda x: rl.check_convert_chr_id(x))

    df_contex['target_key'] =  df_contex['chrom_1st'].astype(str)+ ';' + df_contex['strand_1st'].astype(str)
    df_contex['query_key'] =  df_contex['chrom_2end'].astype(str)+ ';' + df_contex['strand_2end'].astype(str)

    df_contex.to_csv(context_file, index=False)
    return df_contex


def load_occupyed_data(input_occupyed):
    """
    load occupyed data

        Parameters
        ----------
        input_occupyed: input file path

        Returns
        -------
        occupyed_InteLab: Interlab object of occupyed regions
        """
    overlap_handle = open(input_occupyed,'rb')
    occupyed_InteLab = pickle.load(overlap_handle)
    # print(overlap_avg_val)
    overlap_handle.close()
    return occupyed_InteLab


def get_occ_regions(target_key, query_key, target_pos_s,target_pos_e, occupyed_InteLab, query_pos_s,query_pos_e,strand_t, strand_q, block_ends, target_seed_s, target_seed_e, query_seed_s, query_seed_e):
    """
    get_occ_regions

        Parameters
        ----------
        target_key:
        query_key:
        target_pos_s:
        target_pos_e:
        occupyed_InteLab:
        query_pos_s:
        query_pos_e:
        strand_t:
        strand_q:
        block_ends:
        target_seed_s:
        target_seed_e:
        query_seed_s:
        query_seed_e:


        Returns
        -------
        pos_param_t:
        pos_param_q:
        neg_param_t:
        neg_param_q:
        """
    occupied_regions_target = find_occu_overlaps(target_key,
                                                 target_pos_s,target_pos_e,
                                                 occupyed_InteLab)
    occupied_regions_query = find_occu_overlaps(query_key,
                                                query_pos_s,query_pos_e,
                                                occupyed_InteLab)

    # covert occupyed prositons:
    new_occupied_reg_t, pos_occ_list_t = convert_occu_positons(target_pos_s,
                                               target_pos_e,
                                               occupied_regions_target,
                                               strand_t, block_ends,
                                               target_seed_s, target_seed_e)
    new_occupied_reg_q, pos_occ_list_q = convert_occu_positons(query_pos_s,
                                               query_pos_e,
                                               occupied_regions_query,
                                               strand_q, block_ends,
                                               query_seed_s, query_seed_e)

    if pos_occ_list_t:
        pos_param_t = get_neg_pos_intarna_str(pos_occ_list_t, '--tAccConstr=\"b:', target_pos_s, target_pos_e)
    else:
        pos_param_t = ''
    if pos_occ_list_q:
        pos_param_q = get_neg_pos_intarna_str(pos_occ_list_q, '--qAccConstr=\"b:', query_pos_s, query_pos_e)
    else:
        pos_param_q= ''

    neg_param_t = get_neg_pos_intarna_str(new_occupied_reg_t, '--tAccConstr=\"b:', target_pos_s, target_pos_e)
    neg_param_q = get_neg_pos_intarna_str(new_occupied_reg_q, '--qAccConstr=\"b:', query_pos_s, query_pos_e)
    if neg_param_t == '' or neg_param_q == '':
        print('Warining: full sequence is occupyed!')
        full_seq_occ += 1
        return '', '', '', ''
    return pos_param_t, pos_param_q, neg_param_t, neg_param_q


def report(context_not_full,full_seq_occ,no_neg,lost_inst_pos,lost_inst_neg,no_less_sub_opt_pos,no_less_sub_opt_neg):
    """
    Repots number of lost sequences

        Parameters
        ----------
        context_not_full:
        full_seq_occ:
        neg_param:
        lost_inst_pos:
        lost_inst_neg:
        no_less_sub_opt_pos:
        no_less_sub_opt_neg:

        """
    print('####\nContext could not be extended for %i sequences'%context_not_full)
    print('####\nFor %i sequences target and or query is full occupyed'%full_seq_occ)

    print('####\nIntaRNA calls failed:')
    print('%i number of positive IntaRNA calls did not lead to a result'%lost_inst_pos)
    if not no_neg:
        print('%i number of negative IntaRNA calls did not lead to a result'%lost_inst_neg)


    print('####\nNumber of sequences haveing not all suboptimals:')
    print('%i number of positive IntaRNA calls not all suboptimals'%no_less_sub_opt_pos)
    if not no_neg:
        print('%i number of negative IntaRNA calls not all suboptimals'%no_less_sub_opt_neg)


def get_context_file_name(context, pos_occ, block_ends, output_path, experiment_name):
    """
    get_context_file_name

        Parameters
        ----------
        context:
        pos_occ:
        block_ends:
        output_path:
        experiment_name:

        Returns
        -------
        context_file: all no of instaces (trusted rris)
        """
    context_info = '_context_' +  str(context)
    if pos_occ:
        print('%%%%%\nIntaRNA calls including occupied regions for pos and neg instances!!!!\n%%%%%')
        context_info = context_info + '_pos_occ_'
    out_info =  '_block_ends_' +  str(block_ends) + '_'
    context_file = (output_path + experiment_name +  context_info + out_info +
                    'RRI_dataset.csv')
    return context_file, context_info


def get_report_steps(df_contex):
    """
    get_report_steps

        Parameters
        ----------
        df_contex: datafram continig all instances that will be computed

        Returns
        -------
        data_100: all no of instaces (trusted rris)
        data_25: 25 % of number of instances
        data_50: 50 % of number of instances
        date_75: 75 % of number of instances
        """
    data_100 = len(df_contex)
    data_25 = int(data_100*25/100)
    data_50 = int(data_100*50/100)
    date_75 = int(data_100*75/100)
    return data_100, data_25, data_50, date_75

def print_status(index, data_100, data_25, data_50, date_75):
    """
    Report status of the call!

        Parameters
        ----------
        index:
        data_100: all no of instaces (trusted rris)
        data_25: 25 % of number of instances
        data_50: 50 % of number of instances
        date_75: 75 % of number of instances

        """
    if (index+1) == data_100:
        print('***\n full data (%i sequences)\n****' %(data_100))
    elif (index+1) == data_25:
        print('***\n25 percent of the data (%i sequences)\n****' %(data_25))
    elif (index+1) == data_50:
        print('***\n50 percent of the data (%i sequences)\n****' %(data_50))
    elif (index+1) == date_75:
        print('***\n75 percent of the data (%i sequences)\n****' %(date_75))

    # header:
def get_header():
    """
    gets df header and header for IntaRNA call

        Returns
        ----------
        index:
        header: df header
        list_rows_add: output format for IntaRNA call


        """
    list_rows_add = ['score_seq_1st_side', 'score_seq_2end_side',
                     'biotype_region_1st', 'biotype_region_2end',
                     'ID_1st','ID_2end','con_target','con_query',
                     'target_con_s','target_con_e','query_con_s',
                     'query_con_e', 'target_key', 'query_key']
    intaRNA_col_name = 'id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,E_hybrid,ED1,ED2'
    list_intaRNA_col_name = intaRNA_col_name.split(',')
    header = list_intaRNA_col_name + list_rows_add
    return header, list_rows_add


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i1", "--input_rris", action="store", dest="input_rris",
                        required=True,
                        help= "path to file storing all trusted RRIs")
    parser.add_argument("-i2", "--input_occupyed", action="store", dest="input_occupyed",
                        required=True,
                        help= "path to file storing to bindings sides occupying regions")
    parser.add_argument("-d", "--output_path", action="store", dest="output_path",
                        required=True,
                        help= "path output reposetory")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of positve trusted RRIs")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file",
                        required=True, help= "path to 2bit genome file")
    parser.add_argument("-c", "--context",  nargs='?', type=int,
                        dest="context",  default=5,
                        help= "how much context should be added at left an right of the sequence")
    parser.add_argument('--pos_occ', default=False, action='store_true')
    parser.add_argument('--no_pos_occ', dest='pos_occ', action='store_false')
    parser.add_argument("-b", "--block_ends",  nargs='?', type=int,
                        dest="block_ends",  default=0,
                        help= "# nts blocked at the ends of the sequence")
    parser.add_argument("-l", "--chrom_len_file",  action="store", dest="chrom_len_file",
                        required=True,
                        help= "tabular file containing chrom name \t chrom lenght for each chromosome")


    args = parser.parse_args()
    input_rris = args.input_rris
    input_occupyed = args.input_occupyed
    output_path = args.output_path
    experiment_name = args.experiment_name
    genome_file = args.genome_file
    context = args.context
    pos_occ = args.pos_occ
    block_ends = args.block_ends
    chrom_len_file = args.chrom_len_file

    no_sub_opt = 5
    flag_all_neg = True
    flag_all_pos = True
    no_neg = False
    #block_ends = 20
    # pos_occ = False
    #chrom_len_file = '/vol/scratch/data/genomes/mm10.chrom.sizes'
    #chrom_len_file = '/vol/scratch/data/genomes/hg38_Info.tab'

    context_file, context_info = get_context_file_name(context, pos_occ,
                                                       block_ends,
                                                       output_path,
                                                       experiment_name)
    # load occupyed data
    occupyed_InteLab = load_occupyed_data(input_occupyed)


    # Reporting how many instances did not lead to a result
    context_not_full = 0
    lost_inst_pos = 0
    lost_inst_neg = 0
    no_less_sub_opt_pos = 0
    no_less_sub_opt_neg = 0
    full_seq_occ = 0

    ####### Context ###########
    if os.path.isfile(context_file):
        df_contex = pd.read_table(context_file, sep=",")
        print('used existing context file: %s'%(context_file))
    else:
        df_contex = get_context_added(input_rris, output_path, genome_file, context, context_not_full, context_file,chrom_len_file)
        print('***\ncontext is appende pos and negative data generation is starting:\n****')


    ### prepare calling ########################
    data_100, data_25, data_50, date_75 = get_report_steps(df_contex)

    header, list_rows_add = get_header()

    df_pos_data = pd.DataFrame(columns=header)
    df_neg_data = pd.DataFrame(columns=header)

    ### Start producing training instances ###########
    for index, row in df_contex.iterrows():
        # sequences
        target_seq = row['con_target']
        query_seq = row['con_query']

        # potisitons
        target_pos_s = row['target_con_s']
        target_pos_e = row['target_con_e']
        query_pos_s = row['query_con_s']
        query_pos_e = row['query_con_e']

        # postions for seed interaction
        target_seed_s = row['start_1st']
        target_seed_e = row['end_1st']
        query_seed_s = row['start_2end']
        query_seed_e = row['end_2end']

        strand_t = row['strand_1st']
        strand_q = row['strand_2end']

        #t_chr = check_convert_chr_id(row['target_key'].split(';')[0])
        #q_chr = check_convert_chr_id(row['query_key'].split(';')[0])

        #t_id = t_chr + ';' + strand_t
        #q_id = q_chr + ';' + strand_q

        ########### Interacting sequeces
        #print('### RRI Interacting sequeces ########')
        #print('target inter seq: ',row['ineraction_side_1st'])
        #print('target inter seq: ',row['ineraction_side_2end'])


        ########### find places which are occupyed
        #print('### RRI strands ########')
        #print('target side infos: ',strand_t)
        #print('query side infos: ',strand_q)


        ########### find places which are occupyed
        #print('### RRI genome positions ########')
        #print('target side infos: ',target_seed_s,target_seed_e)
        #print('query side infos: ',query_seed_s,query_seed_e)
        #print('### sequence genome positions ########')
        #print('target sequence infos: ',target_pos_s,target_pos_e)
        #print('query sequence infos: ',query_pos_s,query_pos_e)
        #occupied_regions_target = []
        #occupied_regions_query = []

        ##### occupyed regions:

        pos_param_t, pos_param_q, neg_param_t, neg_param_q = get_occ_regions(row['target_key'], row['query_key'], target_pos_s,target_pos_e, occupyed_InteLab, query_pos_s,query_pos_e,strand_t, strand_q, block_ends, target_seed_s, target_seed_e, query_seed_s, query_seed_e)
        if not neg_param_t:
            # this instance will be ignored!
            continue
        #print(occupied_regions_target)

        ####### IntaRNA_call preparation: ######################
        output_columns = ('--outCsvCols id1,start1,end1,id2,start2,end2,' +
                          'subseqDP,hybridDP,E,seedStart1,seedEnd1,' +
                          'seedStart2,seedEnd2,seedE,E_hybrid,ED1,ED2')

        call_general = ('IntaRNA -t ' + target_seq + ' -q ' + query_seq +
                       ' --parameterFile=./IntaRNA_param.txt --outNumber=' +
                       str(no_sub_opt))

        ####POSITIVE DATA##########################
        #### covert occupyed prositons:
        t_seed_pos_new = convert_positions(target_pos_s, target_pos_e,
                                           target_seed_s, target_seed_e,
                                           strand_t)
        q_seed_pos_new = convert_positions(query_pos_s, query_pos_e,
                                           query_seed_s, query_seed_e, strand_q)

        ####pos IntaRNA call
        pos_param = (' --seedQRange='+ str(q_seed_pos_new[0]) + '-' +
                     str(q_seed_pos_new[1]) + ' --seedTRange=' +
                     str(t_seed_pos_new[0]) + '-' +
                     str(t_seed_pos_new[1]) + ' --seedMaxE=0 ')

        if pos_occ:
            pos_param_occ = ' ' + pos_param_t + ' ' + pos_param_q + ' '
            call_pos = call_general + pos_param + pos_param_occ + output_columns
        else:
            call_pos = call_general + pos_param + output_columns

        df_pos_data_old = df_pos_data
        #print('call pos data:\n%s'%call_pos)
        df_pos_data,lost_inst_pos_new, no_less_sub_opt_pos = decode_IntaRNA_call(call_pos,
                                                            lost_inst_pos, row,
                                                            list_rows_add,
                                                            df_pos_data,
                                                            no_sub_opt,
                                                            no_less_sub_opt_pos)

        ####NEGATIV DATA##########################

        if not no_neg:
            neg_param = ' ' + neg_param_t + ' ' + neg_param_q + ' '

            call_neg = call_general + neg_param + output_columns
            #print('call neg data:\n%s'%call_neg)
            #print('lost instaces last call:%i\nthis call:%i'%(lost_inst_pos,lost_inst_pos_new))
            if lost_inst_pos_new == lost_inst_pos:
                df_neg_data, lost_i_neg_new, no_less_sub_opt_neg = decode_IntaRNA_call(call_neg,
                                                             lost_inst_neg, row,
                                                             list_rows_add,
                                                             df_neg_data,
                                                             no_sub_opt,
                                                             no_less_sub_opt_neg)
            elif lost_inst_pos_new > lost_inst_pos and flag_all_neg:
                print('appeded neg data although pos data is not found!')
                print('call pos data:\n%s'%call_pos)
                print('call neg data:\n%s'%call_neg)
                df_neg_data, lost_i_neg_new, no_less_sub_opt_neg = decode_IntaRNA_call(call_neg,
                                                             lost_inst_neg, row,
                                                             list_rows_add,
                                                             df_neg_data,
                                                             no_sub_opt,
                                                             no_less_sub_opt_neg)
                #print('neg instaces:',len(df_neg_data))
                #print('pos instaces:',len(df_pos_data))
            elif lost_inst_pos_new > lost_inst_pos:
                print('negative data not appended because positive not all suboptimals')

            # if negative instance could not be computed check if postive should be ignored
            if not flag_all_pos and lost_i_neg_new > lost_inst_neg:
                df_pos_data = df_pos_data_old
                lost_inst_neg = lost_i_neg_new
            elif lost_i_neg_new > lost_inst_neg:
                lost_inst_neg = lost_i_neg_new
        lost_inst_pos = lost_inst_pos_new


        print_status(index, data_100, data_25, data_50, date_75)

    #print(df_result_neg['start1'], df_result_neg['end1'])
    #print(df_pos_data['start1'],  df_result_neg['end1'])
    result_file = output_path + experiment_name +  context_info
    if no_neg:
        df_pos_data.to_csv(result_file + 'pos.csv', index=False)
    else:
        df_neg_data.to_csv(result_file + 'neg.csv', index=False)
        df_pos_data.to_csv(result_file + 'pos.csv', index=False)

    report(context_not_full,full_seq_occ,no_neg,lost_inst_pos,lost_inst_neg,no_less_sub_opt_pos,no_less_sub_opt_neg)


if __name__ == '__main__':
    main()
