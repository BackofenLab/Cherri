#!/usr/bin/env python
import pandas as pd
import math
import numpy as np
#import matplotlib as mpl
from collections import defaultdict
from interlap import InterLap
import sys
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

def bed_extract_sequences_from_2bit(in_bed, out_fa, in_2bit,
                                    lc_repeats=False,
                                    convert_to_rna=False):
    """
    Extract sequences from genome (provide genome .2bit file).
    twoBitToFa executable needs to be in PATH. Store extracted
    sequences in out_fa.

    !!! COPY !!!!
        Parameters
        ----------

        convert_to_rna:
            If true, read in extracted sequences and convert to RNA.
        lc_repeats:
            If True, do not convert repeat regions to uppercase and output.

        Raises
        ------
        nothing

        Returns
        -------
        seqs_dic
            dictinary holding the sequence ID as key and seqeuce as value

    """
    # Check for twoBitToFa.
    #assert is_tool("twoBitToFa"), "twoBitToFa not in PATH"

    # Run twoBitToFa and check.
    check_cmd = "twoBitToFa"
    if not lc_repeats:
        check_cmd += " -noMask"
    check_cmd += " -bed=" + in_bed + " " + in_2bit + " " + out_fa
    # print(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "twoBitToFa is complaining:\n%s\n%s" %(check_cmd, output)
    #print(output)
    if convert_to_rna:
        # Read in tmp_fa into dictionary (this also converts sequences to RNA).
        seqs_dic = rl.read_fasta_into_dic(out_fa, skip_n_seqs=True)
        # Output RNA sequences.
        #fasta_output_dic(seqs_dic, out_fa, split=True)
    return seqs_dic


def check_context(df, start, end, chrom_end):
    """
    check that the extende contxt is not to short or long!

        Parameters
        ----------
        df: bed df


        Returns
        -------
        df
            df with changed postions

        """
    print('Warning: added context is smaller than 0 for %i instances'%len(df[df.start_1st <=0]))
    df.loc[df.start_1st <= 0, 'start_1st'] = 1

    return df



def get_context(seq_tag, df, out_dir, in_2bit_file, context):
    """
    defining column with ID and empty colums to store the context sequences
    !!! COPY !!!

        Parameters
        ----------
        seq_tag: dataframe
        df: dataframe contining position of the extraction
        out_dir: directory where to store bed and fa file
        in_2bit_file: genome 2bit file
        context: amout of nt that should be added on both sides

        Raises
        ------
        nothing

        Returns
        -------
        df
            colum update dataframe

        """
    out_bed = out_dir + seq_tag + '_out.bed'
    out_fa = out_dir + seq_tag + '_out.fa'
    if seq_tag == 'target':
        df_bed = df[['chrom_1st', 'start_1st', 'end_1st', 'ID1', 'interaction_no', 'strand_1st']].copy()
        #print(df_bed.tail())
        df_bed['chrom_1st'] = df_bed['chrom_1st'].apply(lambda x: rl.check_convert_chr_id(x))
        df_context =  rl.add_context(df_bed, context, 'start_1st', 'end_1st')
        df_context = check_context(df_context, 'start_1st', 'end_1st', 100000000000 )
        # check context!
        col_name = 'con_target'
        col_id = 'ID1'
        df_context_filted = df_context[df_context.chrom_1st != False]
        no_del_entys = len(df_context) - len(df_context_filted)
    elif seq_tag == 'query':
        df_bed = df[['chrom_2end', 'start_2end', 'end_2end', 'ID2', 'interaction_no', 'strand_2end']].copy()
        df_bed['chrom_2end'] = df_bed['chrom_2end'].apply(lambda x: rl.check_convert_chr_id(x))
        df_context =  rl.add_context(df_bed, context, 'start_2end', 'end_2end')
        col_name = 'con_query'
        col_id = 'ID2'
        df_context_filted = df_context[df_context.chrom_2end != False]
        no_del_entys = len(df_context) - len(df_context_filted)
    else:
        print('error: please specify the parameter seq_tag with target or query')
    # delet all 'False' chromosmes of in the df
    print('loost %i instaces because of the Chromosome'%(no_del_entys))
    df_context_filted.to_csv(out_bed, sep="\t", index=False, header=False)
    #df = df_context
    seqs_dic = bed_extract_sequences_from_2bit(out_bed, out_fa, in_2bit_file,lc_repeats=False, convert_to_rna=True)

    for seq_id in seqs_dic:
        #print(seq_id)
        #print(seqs_dic[seq_id])
        df.loc[df[col_id] == seq_id, [col_name]] = seqs_dic[seq_id]

    return df

def find_occu_overlaps(dict_key,s,e,occupied_regions_list, occupyed_InteLab):
    """
    finding all interaction sides for a given sequence by genomic indexes

        Parameters
        ----------
        dict_key: key chrom:strand of the given sequence neede for the occupyed_InteLab
        s: start postion of given sequence
        e: end postion of given sequence
        occupied_regions_list: list contining all occupied positions for the given sequence
        occupyed_InteLab: Interlap object having all interacting side postions

        Raises
        ------
        nothing

        Returns
        -------
        df
            colum update dataframe

        """
    temp_list = list(occupyed_InteLab[dict_key].find((s,e)))
    if temp_list:
        occupied_regions_list.append(temp_list)
    return occupied_regions_list


def convert_positions(s_seq, e_seq, s_side, e_side):
    """
    return postions one based of side in relation to given sequence

        Parameters
        ----------
        s_seq: start postion of given sequence
        e_seq: end postion of given sequence
        s_side: start postion of given occupied side
        e_side: start postion of given occupied side

        Returns
        -------
        new_e_side
            one based end postion of occupied side
        new_s_side
            one based start postion of occupied side
     >>> convert_positions(20,30,22,24)
     3-5
     >>> convert_positions(20,30,12,34)
     1-11

        """
    # check than side are outside the sequence
    if (s_side <= s_seq) or (e_side >= e_seq):
        if (s_side <= s_seq) and (e_side >= e_seq):
            print('Warning: full context is overlaped')
            new_end = e_seq - s_seq + 1
            return  [1, new_end]
        elif s_side <= s_seq:
            s_side = s_seq
        elif e_side >= e_seq:
            e_side = e_seq
    # compute postions
    new_end = e_seq - s_seq + 1
    new_e_side = e_side - s_seq + 1
    new_s_side = s_side - s_seq + 1
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
    join_result_and_infos

        Parameters
        ----------
        df: IntaRNA result df
        lost_inst: list of colum names to add
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
        for col_n in list_rows_add:
            #print(col_n)
            val[col_n] = [row[col_n]]*len(df)

        df_append = pd.DataFrame.from_dict(val)
    #print(df_append)
    df_result = pd.merge(df, df_append, left_index=True, right_index=True)
    return df_result, lost_inst



def convert_occu_positons(seq_s, seq_e, occupied_regions):
    """
    convert the list conining all occupied postions

        Parameters
        ----------
        seq_s: target/query sequence start pos
        seq_e: target/query sequence end pos
        occupied_regions: occupied regions list [(s,e)(s,e)...]

        Returns
        -------
        new_occ_list
            list with updated postions
        """
    new_occ_list = []
    for i in occupied_regions:
        s = i[0]
        e = i[1]
        pos_new = convert_positions(seq_s, seq_e, s, e)
        new_occ_list.append(pos_new)
    return new_occ_list



def decode_IntaRNA_call(call_pos, lost_inst, row, list_rows_add, df_data, no_sub_opt, no_less_sub_opt):
    """
    decode the IntaRNA call

        Parameters
        ----------
        call_pos:
        lost_inst:
        row:
        list_rows_add:
        df_data:
        no_sub_opt:
        no_less_sub_opt:

        Returns
        -------
        new_occ_list
            list with updated postions
        """
    out = rl.call_script(call_pos,reprot_stdout=True)
    df = decode_Intarna_output(out)
    # print(df_pos)
    df = df.reset_index(drop=True)

    df_result_pos, lost_inst_new = join_result_and_infos(df,
                                                         lost_inst,
                                                         row, list_rows_add)
    if lost_inst < lost_inst_new:
        #print('lost instace')
        lost_inst = lost_inst_new
    else:
        #print('instace appended to data')
        # check if found all number of subotpimals
        if no_sub_opt != len(df_result_pos):
            print('Warning IntaRNA could not find %i interaction but %i'%(no_sub_opt, len(df_result_pos)))
            no_less_sub_opt += 1

        df_data = pd.concat([df_data, df_result_pos])
    return df_data, lost_inst


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




    args = parser.parse_args()
    input_rris = args.input_rris
    input_occupyed = args.input_occupyed
    output_path = args.output_path
    experiment_name = args.experiment_name
    genome_file = args.genome_file
    context = args.context
    no_sub_opt = 5


    context_info = '_context_' +  str(context) + '_'
    context_file = (output_path + experiment_name +  context_info +
                    'RRI_dataset.csv')

    # load occupyed data
    overlap_handle = open(input_occupyed,'rb')
    occupyed_InteLab = pickle.load(overlap_handle)
    # print(overlap_avg_val)
    overlap_handle.close()

    context_not_full = 0


    if os.path.isfile(context_file):
        df_contex = pd.read_table(context_file, sep=",")
        #df_pos_RRIs_result = inial_df()
        #df_neg_RRIs_result = inial_df()
        print('used existing context file: %s'%(context_file))
    else:
        df_RRIs = pd.read_table(input_rris, sep=",")

        # adding context by including infors into the df
        df_RRIs = extention_df(df_RRIs)
        df_target = get_context('target', df_RRIs, output_path,
                                genome_file, context)
        #print(df_target)
        df_context_seq = get_context('query', df_target, output_path,
                                     genome_file, context)
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

        df_contex['target_key'] = (''.join(df_contex['chrom_1st'].astype(str)+ ';' + df_contex['strand_1st'].astype(str)))
        df_contex['query_key'] = (''.join(df_contex['chrom_2end'].astype(str)+ ';' + df_contex['strand_2end'].astype(str)))

        #print(df_contex)
        df_contex.to_csv(context_file, index=False)

        print('***\ncontext is appende pos and negative data generation is starting:\n****')



    ## Generate report steps:
    data_100 = len(df_contex)
    data_25 = int(data_100*25/100)
    data_50 = int(data_100*50/100)
    date_75 = int(data_100*75/100)

    lost_inst_pos = 0
    lost_inst_neg = 0
    no_less_sub_opt_pos = 0
    no_less_sub_opt_neg = 0
    # header:
    list_rows_add = ['score_seq_1st_side', 'score_seq_2end_side',
                     'biotype_region_1st', 'biotype_region_2end',
                     'ID_1st','ID_2end','con_target','con_query',
                     'target_con_s','target_con_e','query_con_s',
                     'query_con_e', 'target_key', 'query_key']
    intaRNA_col_name = 'id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,E_hybrid,ED1,ED2'
    list_intaRNA_col_name = intaRNA_col_name.split(',')
    header = list_intaRNA_col_name + list_rows_add

    df_pos_data = pd.DataFrame(columns=header)
    df_neg_data = pd.DataFrame(columns=header)

    for index, row in df_contex.iterrows():
        # sequences
        target_seq = row['con_target']
        query_seq = row['con_query']
        #print(target_seq)
        #print(query_seq)


        # potisitons
        target_pos_s = row['target_con_s']
        target_pos_e = row['target_con_e']
        query_pos_s = row['query_con_s']
        query_pos_e = row['query_con_e']

        #print(target_pos_s)
        #print(target_pos_e)

        # postions for seed interaction
        target_seed_s = row['start_1st']
        target_seed_e = row['end_1st']
        query_seed_s = row['start_2end']
        query_seed_e = row['end_2end']


        ########### find places which are occupyed
        occupied_regions_target = [(target_seed_s,target_seed_e,'rri_side')]
        occupied_regions_query = [(query_seed_s,query_seed_e,'rri_side')]

        occupied_regions_target = find_occu_overlaps(row['target_key'],
                                                     target_pos_s,target_pos_e,
                                                     occupied_regions_target,
                                                     occupyed_InteLab)
        occupied_regions_query = find_occu_overlaps(row['query_key'],
                                                    query_pos_s,query_pos_e,
                                                    occupied_regions_query,
                                                    occupyed_InteLab)

        #print(occupied_regions_target)

        ####### IntaRNA_call preparation: ######################
        output_columns = ('--outCsvCols id1,start1,end1,id2,start2,end2,' +
                          'subseqDP,hybridDP,E,seedStart1,seedEnd1,' +
                          'seedStart2,seedEnd2,seedE,E_hybrid,ED1,ED2')
        # call_general = ('IntaRNA -t ' + target_seq + ' -q ' + query_seq +
        #                ' --outMode C --seedBP 5 --seedMinPu 0 --accW 150' +
        #                ' --acc N --temperature=37 --outMaxE=-5' +
        #                ' --outOverlap=B ')
        call_general = ('IntaRNA -t ' + target_seq + ' -q ' + query_seq +
                       ' --outMode C --seedBP 5 --seedMinPu 0 --accW 150' +
                       ' --acc N --temperature=37 --outMaxE=-5' +
                       ' --outOverlap=B --outNumber=' + str(no_sub_opt) + ' ')

        ####POSITIVE DATA##########################
        #### covert occupyed prositons:
        t_seed_pos_new = convert_positions(target_pos_s, target_pos_e,
                                           target_seed_s, target_seed_e)
        q_seed_pos_new = convert_positions(query_pos_s, query_pos_e,
                                           query_seed_s, query_seed_e)

        ####pos IntaRNA call
        pos_param = (' --seedQRange='+ str(q_seed_pos_new[0]) + '-' +
                     str(q_seed_pos_new[1]) + ' --seedTRange=' +
                     str(t_seed_pos_new[0]) + '-' +
                     str(t_seed_pos_new[1]) + ' ')
        call_pos = call_general + pos_param + output_columns
        # print('call pos data:\n%s'%call_pos)



        df_pos_data, lost_inst_pos = decode_IntaRNA_call(call_pos,
                                                         lost_inst_pos, row,
                                                         list_rows_add,
                                                         df_pos_data,
                                                         no_sub_opt,
                                                         no_less_sub_opt_pos)

        ####NEGATIV DATA##########################
        # covert occupyed prositons:
        new_occupied_reg_t = convert_occu_positons(target_pos_s, target_pos_e,
                                                   occupied_regions_target)
        new_occupied_reg_q = convert_occu_positons(query_pos_s, query_pos_e,
                                                   occupied_regions_query)

        neg_param_t = get_neg_pos_intarna_str(new_occupied_reg_t, ' --tAccConstr=\"b:', target_pos_s, target_pos_e)
        neg_param_q = get_neg_pos_intarna_str(new_occupied_reg_q, ' --qAccConstr=\"b:', query_pos_s, query_pos_e)
        neg_param = ' ' + neg_param_t + ' ' + neg_param_q + ' '

        call_neg = call_general + neg_param + output_columns
        # print('call neg data:\n%s'%call_neg)

        df_neg_data, lost_inst_neg = decode_IntaRNA_call(call_neg,
                                                         lost_inst_neg, row,
                                                         list_rows_add,
                                                         df_neg_data,
                                                         no_sub_opt,
                                                         no_less_sub_opt_neg)


        if (index+1) == data_100:
            print('***\n full data (%i sequences)\n****' %(data_100))
        elif (index+1) == data_25:
            print('***\n25 percent of the data (%i sequences)\n****' %(data_25))
        elif (index+1) == data_50:
            print('***\n50 percent of the data (%i sequences)\n****' %(data_50))
        elif (index+1) == date_75:
            print('***\n75 percent of the data (%i sequences)\n****' %(date_75))

    #print(df_result_neg['start1'], df_result_neg['end1'])
    #print(df_pos_data['start1'],  df_result_neg['end1'])
    result_file = output_path + experiment_name +  context_info
    df_neg_data.to_csv(result_file + '_neg.table', index=False)
    df_pos_data.to_csv(result_file + '_pos.table', index=False)

    #### Report
    print('####\nContext could not be extended for %i sequences'%context_not_full)

    print('####\nIntaRNA calls failed:')
    print('%i number of positive IntaRNA calls did not lead to a result'%lost_inst_pos)
    print('%i number of negative IntaRNA calls did not lead to a result'%lost_inst_neg)


    print('####\nNuber of sequences haveing not all suboptimals:')
    print('%i number of positive IntaRNA calls not all suboptimals'%no_less_sub_opt_pos)
    print('%i number of negative IntaRNA calls not all suboptimals'%no_less_sub_opt_neg)

if __name__ == '__main__':
    main()
