#!/usr/bin/env python
import pandas as pd
import math
import numpy as np
#import matplotlib as mpl
from collections import defaultdict
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


def bed_extract_sequences_from_2bit(in_bed, out_fa, in_2bit,
                                    lc_repeats=False,
                                    convert_to_rna=False):
    """
    Extract sequences from genome (provide genome .2bit file).
    twoBitToFa executable needs to be in PATH. Store extracted
    sequences in out_fa.
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
        seqs_dic = rl.read_fasta_into_dic(out_fa)
        # Output RNA sequences.
        #fasta_output_dic(seqs_dic, out_fa, split=True)
    return seqs_dic


def check_convert_chr_id(chr_id):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

        Parameters
        ----------
        chr_id: chromosme id string

        Raises
        ------
        nothing

        Returns
        -------
        chr_id
            updated chromosme id

    """
    assert chr_id, "given chr_id empty"
    chr_id = str(chr_id)

    if re.search("^chr", chr_id):
        if not re.search("^chr[\dMXY]+$", chr_id):
            chr_id = False
    else:
        # Convert to "chr" IDs.
        if chr_id == "MT":
            chr_id = "M"
        if re.search("^[\dMXY]+$", chr_id):
            chr_id = "chr" + chr_id
        else:
            chr_id = False
    return chr_id


def add_context(df_bed, context, start, end):
    """
    edding the changing the start and end postion of the sequences
    to add context to both sides of the sequences in the dataframe

        Parameters
        ----------
        df_bed: dataframe containing start and end positon
        context: amount of nucleotied
        start: column name of start positons
        end: column name of end positons

        Raises
        ------
        nothing

        Returns
        -------
        df_bed
            datafram with updated postions

        """
    #print(df_bed[start])
    df_bed[start] = df_bed[start] - context
    df_bed[end] = df_bed[end] + context
    #print(df_bed[start])
    return df_bed


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
    df['con_seq_only_target'] = ''
    df['con_seq_only_query'] = ''
    return df

def get_context(seq_tag, df, out_dir, in_2bit_file, context):
    """
    defining column with ID and empty colums to store the context sequences

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
        df_bed['chrom_1st'] = df_bed['chrom_1st'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_1st', 'end_1st')
        col_name = 'con_target'
        col_name_seq_only = 'con_seq_only_target'
        col_id = 'ID1'
    elif seq_tag == 'query':
        df_bed = df[['chrom_2end', 'start_2end', 'end_2end', 'ID2', 'interaction_no', 'strand_2end']].copy()
        df_bed['chrom_2end'] = df_bed['chrom_2end'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_2end', 'end_2end')
        col_name = 'con_query'
        col_name_seq_only = 'con_seq_only_query'
        col_id = 'ID2'
    else:
        print('error: please specify the parameter seq_tag with target or query')
    df_context.to_csv(out_bed, sep="\t", index=False, header=False)
    seqs_dic = bed_extract_sequences_from_2bit(out_bed, out_fa, in_2bit_file,lc_repeats=False, convert_to_rna=True)
    #print()
    for seq_id in seqs_dic:
        #print(seq_id)
        df.loc[df[col_id] == seq_id, [col_name]] = seqs_dic[seq_id]
        context_seq = separate_context_seq(seqs_dic[seq_id], context)
        df.loc[df[col_id] == seq_id, [col_name_seq_only]] = context_seq


    return df

def separate_context_seq(seq, context):
    """
    separate the context from the binding side sequence

        Parameters
        ----------
        seq: sequence
        context: amount of nt added on both sides

        Raises
        ------
        nothing

        Returns
        -------
        context_seq
            context from both sides separated by a '&'

        """
    seq_list = list(seq)
    #print(seq_list)
    con_5prime = seq_list[0:context]
    #print(con_5prime)
    end = len(seq_list)
    end_seq = end - context
    con_3prime = seq_list[end_seq:end]
    #print(con_3prime)
    context_seq = ''.join(con_5prime) + "&" + ''.join(con_3prime)
    #print(context_seq)
    return context_seq

def check_context_extention(df, context, output_path):
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
    df['con_seq_only_target'] = df['con_seq_only_target'].replace('', np.nan)
    df['con_seq_only_query'] = df['con_seq_only_query'].replace('', np.nan)
    #null_col_target = df['con_seq_only_target'].isnull().sum()
    #null_col_query = df['con_seq_only_query'].isnull().sum()

    df_filtered_target = df.dropna(axis=0, subset=['con_seq_only_target'])
    df_filtered_query = df_filtered_target.dropna(axis=0, subset=['con_seq_only_query'])
    # print(len(df))
    # print(len(df_filtered_query))
    null_col = len(df)-len(df_filtered_query)

    if null_col > 0:
        print('warning: context was not fully added for %i RRIs' %(null_col))
    return df_filtered_query






def get_neg_instance(df_pos, df_neg):
    """
    Initializing dataframe with the IntaRNA coulum header

        Parameters
        ----------
        df_pos: dataframe containing the positive RRI InraRNA prediction
        df_neg: dataframe conainging the IntrRNA reusults of the suffled RRI
        sequences of the pos RRI

        Raises
        ------
        nothing

        Returns
        -------
        df_neg_entry
            datafram containing the one neg instance with the closed enery
            to the positive instance

        """
    count_nan_neg = 0
    E_pos_RRI = float(df_pos.E.tolist()[0])
    #print(E_pos_RRI)
    E_neg_RRI_list = [float(i) for i in df_neg.E.tolist()]
    E_neg_RRI_clean_list = [x for x in E_neg_RRI_list if str(x) != 'nan']
    #print(E_neg_RRI_clean_list)


    if len(E_neg_RRI_clean_list) != 0:
        #(lambda x: x*x)(x)
        substrced_list = [abs(x-E_pos_RRI) for x in E_neg_RRI_clean_list]
        #print(substrced_list)
        closest_E = min(E_neg_RRI_clean_list, key=lambda x:abs(x-E_pos_RRI))
    elif len(E_neg_RRI_clean_list) == 0:
        # did not find a negative sample
        closest_E = 'nan'
        count_nan_neg += 1
    #print(closest_E)

    # df_neg = df_neg.astype({'E': 'float'}).dtypes
    df_neg['E'] = df_neg['E'].astype(float)
    #print(df_neg)

    df_neg_entry = df_neg[df_neg['E'] == float(closest_E)]
    #print(df_neg_entry)
    return df_neg_entry, count_nan_neg



def inial_df():
    """
    Initializing dataframe with the IntaRNA coulum header

        Parameters
        ----------


        Raises
        ------
        nothing

        Returns
        -------
        df_inital
            emty df with just the header (columns)

        """
    df_inital = pd.DataFrame(columns=['id1','start1','end1','id2','start2',
                                      'end2','subseqDP','hybridDP','E',
                                      'seedStart1','seedEnd1','seedStart2',
                                      'seedEnd2', 'target', 'query',
                                      'id_target', 'id_query'])
    return df_inital


def Intarna_call(seq1, seq2,df, id_target, id_query):
    """
    Intarna_call

        Parameters
        ----------
        seq1: tartet sequence
        seq2: query sequence
        df: datafame where the IntaRNA result will be appended

        Raises
        ------
        nothing

        Returns
        -------
        df_result
            df incuding IntaRNA result

        """
    #print(df)
    #print(seq1)
    #print(seq2)
    output_columns = '--outCsvCols id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,seedStart1,seedEnd1,seedStart2,seedEnd2'
    call = 'IntaRNA -t ' + seq1 + ' -q ' + seq2 + ' --outMode C --seedBP 5 --seedMinPu 0 --accW 150 --acc N --temperature=37 ' + output_columns
    # print(call)

    process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
    #print(call)
    for idx, line in enumerate(process.stdout):
        line = line.decode("utf-8").strip().split(';')
        #line = line.strip()
        #line = line.split(';')
        #print(line)
        if idx == 0:
            col = line
            result = 'empty'
        elif idx == 1:
            values = line
            result = 'exists'
        else:
            print('error IntaRNA output has unexpected number of rows')
        #print(line)
    #print(result)
    if result == 'empty':
        no_values = ['nan','nan','nan','nan','nan','nan','nan',
                     'nan','nan','nan','nan','nan','nan']
        #print(col)
        df_one_result = pd.DataFrame([no_values], columns=col)
    elif result == 'exists':
        df_one_result = pd.DataFrame([values], columns=col)

    df_one_result['target'] = seq1
    df_one_result['query'] = seq2
    df_one_result['id_target'] = id_target
    df_one_result['id_query'] = id_query
    df_result = pd.concat([df, df_one_result])

    #time.sleep(15)
    #assert os.path.isfile(temp_out_cvs_path), 'temp_out_cvs_path not found'
    #if os.stat(temp_out_cvs_path).st_size == 0:
        #print('file is empty')
        #print(temp_out_cvs_path)
        #sys.exit()

    #f = open(temp_out_cvs_path, "r")
    #print(f.read())
    #df_one_result = pd.read_table(temp_out_cvs_path, sep=";")
    #df_result = pd.concat([df, df_one_result])
    #os.remove(temp_out_cvs_path)
    #else:
        #df_result = df
    return df_result



def predict_hybrid_for_neg_seq(shuffled_target_list, shuffled_query_list, id_target, id_query):
    """
    predict hybrid for neg seq

        Parameters
        ----------
        shuffled_target_list: list of negative sequences

        Raises
        ------
        nothing

        Returns
        -------
        df_neg
            df having IntRNA prediction results

        """
    df_initial_neg_result = inial_df()
    for idx2, neg_target in enumerate(shuffled_target_list):
        neg_query = shuffled_query_list[idx2]
        if idx2 == 0:
            #print('if')
            df_neg = Intarna_call(neg_target, neg_query, df_initial_neg_result, id_target, id_query)
        else:
            #print('else')
            df_neg = Intarna_call(neg_target, neg_query, df_neg, id_target, id_query)

    return df_neg


def choose_shuffling(hybrid_seq, hybrid, target, query,
                               shuffle_no_seq, kind_of_shuffel):
    """
    choose shuffling method. This sequence is calling the shuffleing method
    based on the kind of shuffling parameter.

        Parameters
        ----------
        shuffled_target_list: list of negative sequences

        Raises
        ------
        nothing

        Returns
        -------
        shuffled_target_list
            shuffled target sequence list
        shuffled_query_list
            shuffled target sequence list
        """
    if kind_of_shuffel == '3':
        shuffled_target_list, shuffled_query_list = rl.bp_suffeling(hybrid_seq,
                                                         hybrid, shuffle_no_seq)
    elif kind_of_shuffel == '2' or kind_of_shuffel == '1':
        shuffled_target_list = rl.shuffle_sequence(target, shuffle_no_seq,
                                                kind_of_shuffel)
        shuffled_query_list = rl.shuffle_sequence(query, shuffle_no_seq,
                                               kind_of_shuffel)
    else:
        print('Error: please provied a kind of shuffleing 1, 2 or 3')

    return shuffled_target_list, shuffled_query_list

def get_shuffled_context(context, kind_of_shuffel):
    """
    get_shuffled_context

        Parameters
        ----------
        context: context sequenes sepatated with a &
        kind_of_shuffel: suffeling kind [1 or 2 or 4]
            4 means no shuffling

        Raises
        ------
        nothing

        Returns
        -------
        con_shuffl_1
            shuffled 5' context sequence string
        con_shuffl_2
            shuffled 3' context sequence string
        """

    context_list = context.split("&")
    if kind_of_shuffel == 4:
        con_shuffl_1 = context_list[0]
        con_shuffl_2 = context_list[1]
    else:
        con_shuffl_1_list = rl.shuffle_sequence(context_list[0], 1, kind_of_shuffel)
        con_shuffl_2_list = rl.shuffle_sequence(context_list[1], 1, kind_of_shuffel)
        con_shuffl_1 = con_shuffl_1_list[0]
        con_shuffl_2 = con_shuffl_2_list[0]
    return con_shuffl_1, con_shuffl_2


def get_seq_with_context(shuffled_list, con_sh_1, con_sh_2):
    """
    get_shuffled_context

        Parameters
        ----------
        shuffled_list: list of shuffled sequences
        con_sh_1: suffeling kind [1 or 2]

        Raises
        ------
        nothing

        Returns
        -------
        seq_con_list
            shuffled sequences including sequencing
        """
    seq_con_list = []
    for seq in shuffled_list:
        seq_con = con_sh_1 + seq + con_sh_2
        seq_con_list.append(seq_con)

    return seq_con_list



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                        required=True,
                        help= "path to file storing all positve trusted RRIs")
    parser.add_argument("-d", "--output_path", action="store", dest="output_path",
                        required=True,
                        help= "path output reposetory")
    parser.add_argument("-n", "--experiment_name", action="store",
                        dest="experiment_name", required=True,
                        help= "name of the datasoruce of positve trusted RRIs")
    parser.add_argument("-k", "--kind_of_shuffel",  nargs='?',
                        dest="kind_of_shuffel",  default=2,
                        help= "seqence mononucleotide (1) or sequence denucleotide (2) or bp mononucleotide (3) shuffling")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file",
                        required=True, help= "path to 2bit genome file")
    parser.add_argument("-s", "--shuffle_no_seq",  nargs='?',
                        dest="shuffle_no_seq",  default=5, type=int,
                        help= "how often is the positive sequence shuffled")
    parser.add_argument("-cm", "--context_method",  nargs='?',
                        dest="context_method",  default='non',
                        help= "select the context method  if context should not be added (non), if it should be shuffled sepatatly (separat), or together (together) with the sequence")
    parser.add_argument("-c", "--context",  nargs='?', type=int,
                        dest="context",  default=5,
                        help= "how much context should be added at left an right of the sequence")




    args = parser.parse_args()
    input_file = args.input_file
    output_path = args.output_path
    experiment_name = args.experiment_name
    kind_of_shuffel = args.kind_of_shuffel
    genome_file = args.genome_file
    shuffle_no_seq = args.shuffle_no_seq
    context_method = args.context_method
    context = args.context
    #output_path = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/'
    #print(type(kind_of_shuffel))
    # context_method = 'non', 'separat', 'together'
    #context_method = 'together'
    #context = 4
    #out_dir= '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/output/'
    #in_2bit_file = '/home/teresa/Dokumente/RNA_RNA_interaction_evaluation/data/genomes/hg38_UCSC_20210318.2bit'

    #print(type(shuffle_no_seq))


    context_info = '_context_method_' + context_method + '_shuffling_method_' + kind_of_shuffel + '_with_' + str(context) + '_context_'
    context_file = output_path + experiment_name +  context_info + 'RRI_dataset.csv'

    print('context method %s with shuffling method %s' %(context_method, kind_of_shuffel))

    if os.path.isfile(context_file):
        df_filted_RRIs = pd.read_table(context_file, sep=",")
        df_pos_RRIs_result = inial_df()
        df_neg_RRIs_result = inial_df()
        print('used existing context file: %s'%(context_file))
    else:
        df_pos_RRIs_result = inial_df()
        df_neg_RRIs_result = inial_df()

        df_RRIs = pd.read_table(input_file, sep=",")

        # adding context by including infors into the df
        df_RRIs = extention_df(df_RRIs)
        df_target = get_context('target', df_RRIs, output_path, genome_file, context)
        df_context = get_context('query', df_target, output_path, genome_file, context)

        df_filted_RRIs = check_context_extention(df_context, context, output_path)
        #### context added df saved!
        df_filted_RRIs.to_csv(context_file, index=False)

        print('***\ncontext is appende negative data generation by suffeling is starting:\n****')


    if context_method == 'separat':
        if int(kind_of_shuffel) != 1 and int(kind_of_shuffel) != 2:
            print('error: for shuffling method separat please only choose 1 or 2 as kind of shuffling')
            sys.exit()

    # print('Kind of shuffeling: %i' % (int(kind_of_shuffel)))


    ### Generate report steps:
    data_100 = len(df_filted_RRIs)
    data_25 = int(data_100*25/100)
    data_50 = int(data_100*50/100)
    date_75 = int(data_100*75/100)



    for index, row in df_filted_RRIs.iterrows():
        target = row['ineraction_side_1st']
        query = row['ineraction_side_2end']
        hybrid_seq = target + '&' + query
        hybrid = row['IntaRNA_prediction']

        #print(row)

        if (index+1) == data_100:
            print('***\n full data (%i sequences)\n****' %(data_100))
        elif (index+1) == data_25:
            print('***\n25 percent of the data (%i sequences)\n****' %(data_25))
        elif (index+1) == data_50:
            print('***\n50 percent of the data (%i sequences)\n****' %(data_50))
        elif (index+1) == date_75:
            print('***\n75 percent of the data (%i sequences)\n****' %(date_75))

        #pos sequence
        df_initial_pos_result = inial_df()
        df_pos = Intarna_call(target, query, df_initial_pos_result, row['ID1'], row['ID2'])



        if context_method == 'non':
            shuffled_target_list, shuffled_query_list = choose_shuffling(hybrid_seq, hybrid,
                                                                         target, query,
                                                                         shuffle_no_seq,
                                                                         kind_of_shuffel)
        elif context_method == 'separat':
            seq_sh_target_list, seq_sh_query_list = choose_shuffling(hybrid_seq, hybrid,
                                                                         target, query,
                                                                         shuffle_no_seq,
                                                                         kind_of_shuffel)
            #print(row['con_seq_only_target'])
            target_con_sh_1, target_con_sh_2 = get_shuffled_context(row['con_seq_only_target'], kind_of_shuffel)
            query_con_sh_1, query_con_sh_2 = get_shuffled_context(row['con_seq_only_query'], kind_of_shuffel)
            # concatinate context to sequences:
            #print('shufled con 1: %s and 2: %s'%(target_con_sh_1, target_con_sh_2))
            shuffled_target_list = get_seq_with_context(seq_sh_target_list, target_con_sh_1, target_con_sh_2)
            #print(shuffled_target_list[0])
            shuffled_query_list = get_seq_with_context(seq_sh_query_list, query_con_sh_1, query_con_sh_2)

        elif context_method == 'together':
            # get sequences including context and than apply
            #print(row['con_seq_only_target'])
            target_con_sh_1, target_con_sh_2 = get_shuffled_context(row['con_seq_only_target'], 4)
            query_con_sh_1, query_con_sh_2 = get_shuffled_context(row['con_seq_only_query'], 4)
            #print('not shufled con 1: %s and 2: %s'%(target_con_sh_1, target_con_sh_2))

            target_con = target_con_sh_1 + target + target_con_sh_2
            query_con = query_con_sh_1 + query + query_con_sh_2
            hybrid_seq_con = target_con + '&' + query_con
            #print(target_con)
            #print(query_con)
            #print(hybrid_seq_con)
            hybrid_target, hybrid_query = get_shuffled_context(hybrid, 4)
            con_hybrid = "." * context
            hybrid_con = con_hybrid + hybrid_target + con_hybrid + '&' +con_hybrid + hybrid_query + con_hybrid
            #print(hybrid_con)
            shuffled_target_list, shuffled_query_list = choose_shuffling(hybrid_seq_con, hybrid_con,
                                                                         target_con, query_con,
                                                                         shuffle_no_seq,
                                                                         kind_of_shuffel)
        else:
            print('error: please specify only non, separat or together as context method')


##############Call IntaRNA to select the negevie sequence#######################
        #print(df_pos)
        df_neg =  predict_hybrid_for_neg_seq(shuffled_target_list, shuffled_query_list, row['ID1'], row['ID2'])
        #print(df_neg)

#########select the negativ instance closes to the pos energy###################
        df_neg_entry, count_nan_neg = get_neg_instance(df_pos, df_neg)
        #print(df_neg_entry)



        #### save positive and negativ instance in result df
        df_pos_RRIs_result = pd.concat([df_pos_RRIs_result, df_pos])
        df_neg_RRIs_result = pd.concat([df_neg_RRIs_result, df_neg_entry])

    ################################################################

# add all clums of the input data
    #print('negative entry:')
    #print(df_neg_RRIs_result.info())
    #result = pd.concat([df_neg_entry, row], axis=1, keys=['df_neg_entry', 'row'], names=['id_target', 'ID1'])
    df_neg_RRIs_all_result = pd.merge(df_neg_RRIs_result, df_filted_RRIs,  how='left', left_on=['id_target','id_query'], right_on = ['ID1','ID2'])
    df_pos_RRIs_all_result = pd.merge(df_pos_RRIs_result, df_filted_RRIs,  how='left', left_on=['id_target','id_query'], right_on = ['ID1','ID2'])

    #df_filted_RRIs.loc[index]
    #print('appended data:')
    #print(df_neg_RRIs_all_result.info())



    df_neg_RRIs_all_result.to_csv(output_path + experiment_name +  context_info + '_neg_RRI_dataset.csv', index=False)
    df_pos_RRIs_all_result.to_csv(output_path + experiment_name + context_info +'_pos_RRI_dataset.csv', index=False)

    if count_nan_neg > 0:
        print('for %i postivie instances no negative instance was found' % count_nan_neg)


if __name__ == '__main__':
    main()
