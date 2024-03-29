#!/usr/bin/env python3

#general imports

import pandas as pd
import subprocess
import random
import re
import os
from interlap import InterLap
from collections import defaultdict
# training imports
import csv
import numpy as np
#import pandas_profiling
import sklearn as sk
from sklearn import model_selection
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import (RandomForestClassifier)
from sklearn.linear_model import (LogisticRegression)
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import (KNeighborsClassifier)
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import f1_score

#from sklearn.linear_model import LogisticRegression
#import xgboost
import pickle
from sklearn.linear_model import Lasso
from scipy.sparse import csr_matrix, vstack, hstack, load_npz, save_npz
partner =  {a:b for a,b in zip("({[<",")}]>")  }
import eden.graph as eg
import networkx as nx
from ubergauss import tools
from ubergauss.tools import loadfile, dumpfile
from lmz import *
from itertools import compress
import biofilm.algo.feature_selection as featsel
import wget
from Bio import SeqIO
import json

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ CheRRI library ~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Author: muellert [at] informatik.uni-freiburg.de

~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

python3 -m doctest lib.py

"""


################################################################################

def read_chira_data(in_file, header='no', separater="\t"):
    """
    Read RRI tabular file and convert to a dataframe including a header

        Parameters
        ----------
        in_file : tabular file with output of ChiRA RRI results

        Raises
        ------
        nothing

        Returns
        -------
        df_interactions
            dataframe listing all interactions

        """

    # inclued header
    if header == 'no':

        header_line = ['#reads','chrom_1st','start_1st','end_1st', 'strand_1st',
                'chrom_2end','start_2end','end_2end', 'strand_2end',
                'interaction_site_1st', 'interaction_site_2end',
                'IntaRNA_prediction', 'energy',
                'seq_1st_interaction_site', 'seq_2end_interaction_site',
                'start_interaction',
                'chrom_seq_1st_site', 'start_seq_1st_site',
                'stop_seq_1st_site','strand_seq_1st_site',
                'chrom_seq_2end_site', 'start_seq_2end_site',
                'stop_seq_2end_site','strand_seq_2end_site',
                'TPM_seq_1st_site', 'TPM_seq_2end_site', 'TPM_summary',
                'score_seq_1st_site', 'score_seq_2end_site','score_product',
                'biotype_region_1st', 'biotype_region_2end', 'ID_1st','ID_2end']
        df_interactions = pd.read_table(in_file, header=None, sep=separater,
                                low_memory=False, names=header_line)
        #df_interactions = pd.read_table(in_file, header=0, sep=separater)
        #df_temp.columns = header_line
        #df_interactions = df_temp
        #print(df_temp.to_numpy)
        #print(len(header_line))
        #print(df_temp)
        #print(df_interactions.info())

        # df_interactions = pd.DataFrame(df_temp.to_numpy, columns=header_line)
    elif header == 'yes':
        df_interactions = pd.read_table(in_file, sep=separater,
                                        low_memory=False)
    return df_interactions


################################################################################

def get_input_positions_df(df_interactions, needed_header):
    """
    Create a dataframe that contains the interaction information for the user provided format

        Parameters
        ----------
        df_interactions : user defined input in a dataframe
        needed_header: header line required by a user defined input setting


        Returns
        -------
        df_pos
            dataframe containgin the columns matching the needed_header

        """
    if len(df_interactions.columns) == 8:
        assert set(list(df_interactions.columns)) == set(needed_header), 'You input header line is not fitting the needed format'
        df_pos = df_interactions
    elif len(df_interactions.columns) > 8 and len(df_interactions.columns) < 34:
        df_pos = df_interactions[needed_header]
        # print(df_pos)
    else:
        assert len(df_interactions.columns) >= 34, 'Unexpected number of input columns'
    return df_pos

################################################################################

def build_df_and_check_input(df_in, header,NA_list):

    """
    Build the a dataframe having all header columns either with NA list or if existing data stored in df_in

        Parameters
        ----------
        df_in : user defined input in a dataframe
        header: header line of the df that should be build


        Returns
        -------
        df_pos
            dataframe with columns matching the needed_header

        """
    if set(header).issubset(set(list(df_in.columns))):
        df_out = df_interactions[header]
    else:
        data_dict = {}
        for col in header:
            if col in list(df_in.columns):
                data_dict[col] = df_in[col]
            else:
                data_dict[col] = NA_list
        df_out = pd.DataFrame(data_dict)

    return df_out

################################################################################

def set_score_if_not_defined(df_in, df_score, col, score_val=1.0):

    """
    Set a score_val if no score was defined already by the user

        Parameters
        ----------
        df_in : user defined input in a dataframe
        df_score: dataframe containing score information
        col: column name that should be checked and maybe changed
        score_val: score value which should be higher than score_th


        Returns
        -------
        df_score
            dataframe containing changed scores if not defined by the user

    >>> df_in_1 = pd.DataFrame({'score':[0.3,0.3,1,1],'b':[1,1,1,1]})
    >>> df_score_1 = pd.DataFrame({'score':[0.3,0.3,1,1],'b':[1,1,1,1]})
    >>> set_score_if_not_defined(df_in_1, df_score_1, 'score')
       score  b
    0    0.3  1
    1    0.3  1
    2    1.0  1
    3    1.0  1

    >>> df_in_2 = pd.DataFrame({'a':['a','b','c','d'],'b':[1,1,1,1]})
    >>> df_score_2 = pd.DataFrame({'score':['NA','NA','NA','NA'],'b':[1,1,1,1]})
    >>> set_score_if_not_defined(df_in_2, df_score_2, 'score')
       score  b
    0    1.0  1
    1    1.0  1
    2    1.0  1
    3    1.0  1
        """
    if not col in list(df_in.columns):
        df_score[col] = [score_val]*len(df_in)

    return df_score



################################################################################

def check_file_type(in_file):
    """
    Checks file type of the input file to find the delimiter

        Parameters
        ----------
        in_file : inupt file



        Returns
        -------
        delimiter
            type of delimiter

        """
    delimiter = ''
    if in_file.endswith('.tabular'):
        delimiter="\t"
    elif in_file.endswith('.table'):
        delimiter="\t"
    elif in_file.endswith('.csv'):
        delimiter = ","
    else:
        assert delimiter != '', 'Unknown file ending. Use .csv or .tabular'
    return delimiter

################################################################################



def check_input_data(in_file, delimiter):
    """
    Function checks if input data originate form ChiRA or is
    the user defined data.
    If it is a user defined input data is converted to ChiRA format

        Parameters
        ----------
        in_file : tabular file with output of C RRI results
        delimiter : ',' or '\t'


        Returns
        -------
        df_user_in
            user defined input data

        """

    # inclued header
    header_line = ['#reads','chrom_1st','start_1st','end_1st', 'strand_1st',
                  'chrom_2end','start_2end','end_2end', 'strand_2end',
                  'interaction_site_1st', 'interaction_site_2end',
                  'IntaRNA_prediction', 'energy',
                  'seq_1st_interaction_site', 'seq_2end_interaction_site',
                  'start_interaction', 'chrom_seq_1st_site',
                  'start_seq_1st_site', 'stop_seq_1st_site',
                  'strand_seq_1st_site', 'chrom_seq_2end_site',
                  'start_seq_2end_site', 'stop_seq_2end_site',
                  'strand_seq_2end_site', 'TPM_seq_1st_site',
                  'TPM_seq_2end_site', 'TPM_summary', 'score_seq_1st_site',
                  'score_seq_2end_site','score_product', 'biotype_region_1st',
                  'biotype_region_2end', 'ID_1st','ID_2end']

    needed_header = ['chrom_1st','start_1st','end_1st','strand_1st',
                    'chrom_2end','start_2end','end_2end','strand_2end']

    prediction_info_h = ['interaction_site_1st', 'interaction_site_2end',
                        'IntaRNA_prediction', 'energy',
                        'seq_1st_interaction_site',
                        'seq_2end_interaction_site', 'start_interaction']
    prediciton_pos_h = ['chrom_seq_1st_site', 'start_seq_1st_site',
                       'stop_seq_1st_site','strand_seq_1st_site',
                       'chrom_seq_2end_site', 'start_seq_2end_site',
                       'stop_seq_2end_site','strand_seq_2end_site']

    score_h = ['TPM_seq_1st_site', 'TPM_seq_2end_site', 'TPM_summary',
              'score_seq_1st_site', 'score_seq_2end_site','score_product',
              'biotype_region_1st', 'biotype_region_2end', 'ID_1st',
              'ID_2end']


    # check input format
    df_interactions = pd.read_table(in_file, sep=delimiter, low_memory=False)

    #print(df_interactions.info())
    num_interaction = len(df_interactions)
    #print(num_interaction)
    NA_list = ['NA']*num_interaction

    df_reads = pd.DataFrame([1]*num_interaction, columns=['#reads'])


    # check if length is shorter than 34
    if len(df_interactions.columns) == 34 and 'chrom_1st' not in list(df_interactions.columns):
        print('Detected input data in format of ChiRA without header')
        df_interactions = pd.read_table(in_file, header=None, sep=delimiter,
                                low_memory=False, names=header_line)
        return df_interactions
    else:
        print('Convert user input into ChiRA input format')
        df_pos = get_input_positions_df(df_interactions, needed_header)


    #build intaRNA_interaction

    df_pred_i = build_df_and_check_input(df_interactions,
                                            prediction_info_h, NA_list)
    #print(df_pred_i)
    #print(df_pos)
    if not 'chrom_seq_1st_site' in set(df_interactions.columns):
        df_pred_pos = df_pos.copy()
        df_pred_pos.columns = prediciton_pos_h
    else:
        dict_pred_pos = {}
        for col in prediciton_pos_h:
            assert col in df_interactions.columns, 'please provide all predicted interaction positions'
            dict_pred_pos[col] = df_interactions[col]
        df_pred_pos = pd.DataFrame(data=dict_pred_pos)
    #print(df_pos)
    # print(df_pred_pos)

    #df_pred_pos = build_df_and_check_input(df_interactions,
    #                                             prediciton_pos_h, NA_list)
    #print(df_pred_pos)

    df_score_temp = build_df_and_check_input(df_interactions, score_h, NA_list)

    df_score = set_score_if_not_defined(df_interactions, df_score_temp,
                                        'score_seq_1st_site')
    df_score = set_score_if_not_defined(df_interactions, df_score_temp,
                                        'score_seq_2end_site')

    #print(df_reads.info())
    #print(df_pos.info())
    #print(df_pred_i.info())
    #print(df_pred_pos.info())
    #print(df_score.info())

    df_user_in = pd.concat([df_reads, df_pos, df_pred_i, df_pred_pos, df_score],
                           axis=1)
    #print(df_user_in)

    assert len(df_user_in.columns) == 34, 'Wrongly converted user input'

    return df_user_in


################################################################################

def filter_score(df, score_th):
    """
    Filter dataframe for instances with a score of 1

        Parameters
        ----------
        df : df including the containing all RRIs (Interactions)
        score_th: threshold for the expectation maximization score of ChiRA


        Returns
        -------
        df_filterd
            dataframes filter for a score smaller equal threshold of both
            interacting partners

    >>> data = {'score_seq_1st_site':[0.3, 1, 0.9, 0.4],
    ...         'score_seq_2end_site':[0.7, 1, 0.2, 0.5]}
    >>> df = pd.DataFrame(data)
    >>> filter_score(df, 1)
       score_seq_1st_site  score_seq_2end_site
    1                 1.0                  1.0

            """
    # filter input for score_seq_1st_site and score_seq_2end_site == 1
    df_filterd = df[(df.score_seq_1st_site >= score_th) &
                    (df.score_seq_2end_site >= score_th)]
    #df_interactions_single_mapped
    return df_filterd

################################################################################

def delet_empty_col(df):
    """
    Filters out all rows with a empty column entry in a dataframe

        Parameters
        ----------
        df : df including the containing all RRIs
        col_name: name of to be filtered column


        Returns
        -------
        df_filtered
            dataframes without empty column entries


    >>> df = pd.DataFrame({'a':[0.3,'',1,''],'b':[1,1,1,1]})
    >>> delet_empty_col(df)
         a  b
    0  0.3  1
    2  1.0  1
            """
    df_temp = df.replace('', np.nan)
    df_filtered = df_temp.dropna()
    return df_filtered

################################################################################

def call_script(call,reprot_stdout=False,asset_err=True):
    """
    Starts a subprocess to call a script and checks for errors.


        Parameters
        ----------
        call : cmd command

        Returns
        -------
        out
            if reprot_stdout set True returns the stdout (default: False)

    """
    process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
    #process.wait()
    out, err = process.communicate()
    #print(err.decode('utf-8'))
    error = err.decode('utf-8')
    if asset_err:
        assert not error, "script is complaining:\n%s\n%s" %(call, error)
    else:
        if error != "":
            print(f'ERROR or WARNING:\n {error}\nFor call: {call}')
    if reprot_stdout == True:
        return out
        #for line in out.readlines():
            #print(line)


################################################################################

def calculate_overlap(s1,e1,s2,e2,len_flag=False):
    """
    Building for each replicate a inter object

        Parameters
        ----------
        s1: start of one sequence of the first replicate
        e1: end of one sequence of the first replicate
        s2: start of one sequence of the current replicate
        e2: end of one sequence of the current replicate

        Returns
        -------
        combined_overlap
            the combined overlap of sequence 1 and sequence 2

        """
    # print(s1,e1,s2,e2)
    if s1 <= s2:
        s_overlap = s2
        if e1 <= e2:
            e_overlap = e1
        elif e2 < e1:
            e_overlap = e2
        else:
            print('Error: something is not overlapping hire')
    elif s2 < s1:
        s_overlap = s1
        if e1 <= e2:
            e_overlap = e1
        elif e2 < e1:
            e_overlap = e2
        else:
            print('Error: something is not overlapping hire')
    overlap_len = e_overlap - s_overlap +1
    seq1_len = e1 - s1 + 1
    seq2_len = e2 - s2 + 1
    overlap_seq1 = overlap_len/seq1_len
    # print(overlap_seq1)
    overlap_seq2 = overlap_len/seq2_len
    # print(overlap_seq2)
    # combined_overlap = (overlap_seq1 + overlap_seq2)/2
    # select overlap of shorter sequence:
    combined_overlap = max([overlap_seq1, overlap_seq2])
    # print(combined_overlap)
    if len_flag:
        return overlap_len
    else:
        return combined_overlap


################################################################################

def get_chrom_list_no_numbers(df_interactions, chrom):
    """
    Generates a unique list of chromosomes or conticts for both interaction
    partners

        Parameters
        ----------
        df_interactions : df including the filtered RRIs
        chrom : string indicating from which seq the chromosome is


        Returns
        -------
        sort_list_chrom
            sorted list of unique chromosomes or contics which are not a number
            and present in the input data frame

        """
    chrom_list = df_interactions[chrom].unique().tolist()
    #convert all values to string in case it is not
    new_list = [str(el) for idx,el in enumerate(chrom_list)]
    sort_list_chrom = sorted(new_list)

    return sort_list_chrom


################################################################################

def get_list_chrom(df):
    """
    Generates a unique list of chromosomes or contics for both interaction
    partners

        Parameters
        ----------
        df : dataframe including the filtered RRIs


        Returns
        -------
        sort_list_chrom
            sorted list of unique chromosomes or contics which are not a number
            and present in the input data frame

        """
    chrom1_list = get_chrom_list_no_numbers(df, 'chrom_1st')
    chrom2_list = get_chrom_list_no_numbers(df, 'chrom_2end')
    list_chrom_no_int = list(set().union(chrom1_list,chrom2_list))
    sort_list_chrom = sorted(list_chrom_no_int)
    return sort_list_chrom


################################################################################
### functions context

def check_convert_chr_id(chr_id):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertible.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

        Parameters
        ----------
        chr_id: chromosme id string


        Returns
        -------
        chr_id
            updated chromosome id

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


################################################################################

def add_context(df_bed, context, start, end):
    """
    edding the changing the start and end position of the sequences
    to add context to both sites of the sequences in the dataframe

        Parameters
        ----------
        df_bed: dataframe containing start and end position
        context: amount of nucleotide
        start: column name of start position
        end: column name of end position


        Returns
        -------
        df_bed
            dataframe with updated positions

        """
    #print(df_bed[start])
    df_bed[start] = df_bed[start] - context
    df_bed[end] = df_bed[end] + context

    #print(df_bed[start])
    return df_bed


################################################################################

def bed_extract_sequences_from_fasta(in_bed, out_fa, in_genome_fasta,report=1,
                                     all_uc=False,
                                     skip_data_id="set",
                                     skip_n_seqs=True):
    """
    Extract sequences from genome (provide genome FASTA file).
    fastaFromBed executable needs to be in PATH. Store extracted
    sequences in out_fa.

        Parameters
        ----------
        in_bed:
            Input BED file
        out_fa:
            Path to output FASTA file
        in_genome_fasta:
            Genome in FASTA format
        report:
            report the sequence containing N's (1) or
            how many sequences are skipped because of a N (2)


        Returns
        -------
        seqs_dic
            dictionary holding the sequence ID as key and sequence as value

    """
    # Run fastaFromBed and check.
    getfasta_cmd = " ".join(['fastaFromBed', '-s', '-nameOnly',
                        '-fi', in_genome_fasta,
                        '-bed', in_bed,
                        '-fo', out_fa])
    output = subprocess.getoutput(getfasta_cmd)
    # print('fastaFromBed:' + output)
    error = False
    #if output:
    #    error = True
    #assert error == False, "fastaFromBed is complaining:\n%s\n%s" %(getfasta_cmd, output)

    seqs_dic = defaultdict()
    fa_seq = SeqIO.parse(open(out_fa), 'fasta')
    for record in fa_seq:
        # ids of the form chr1:1602:1618:+(+). Strip last 3 characters
        seqs_dic[record.id[:-3]] = str(record.seq).upper().replace('T', 'U')

    # If sequences with N nucleotide should be skipped.
    c_skipped_n_ids = 0
    del_ids = []
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        if re.search("N", seq, re.I):
            if report == 1:
                print (f'WARNING: sequence with seq_id {seq_id} in file {out_fa} contains N nucleotides. Discarding sequence ... ')
            c_skipped_n_ids += 1
            del_ids.append(seq_id)
    for seq_id in del_ids:
        del seqs_dic[seq_id]
    assert seqs_dic, "no sequences remaining after deleting N containing sequences (input FASTA file)"
    if c_skipped_n_ids:
        if report == 2:
            print("# of N-containing %s regions discarded:  %i" %(skip_data_id, c_skipped_n_ids))

    # After converting to RNA and filtering 'N' sequences write again to the FASTA files
    with open(out_fa, "w") as fh_outfa:
        for seq_id, seq in seqs_dic.items():
            fh_outfa.write(">" + seq_id + "\n" + seq + "\n")
    return seqs_dic


################################################################################

def check_context(df, seq_tag, chrom_dict):
    """
    check that the extended context is not to short or long!

        Parameters
        ----------
        df: df with added context
        seq_tag: target or query
        chrom_dict: dictionary storing chromosome length {name->len,...}


        Returns
        -------
        df
            df with changed positions
        count
            countes the sequences where exculded due to false chrom name

    >>> chrom_dict = {'chr1':60,
    ...         'chr2':80}
    >>> datat = {'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'chr1', 'False', 'chr1'],
    ...         'end_1st':[30, 60, 70, 70]}
    >>> dataq = {'start_2end':[0, -40, 40, -2],
    ...         'end_2end':[30, 60, 70, 70],
    ...         'chrom_2end':['chr1', 'chr1', 'chr1', 'chr1']}
    >>> dft = pd.DataFrame(datat)
    >>> dfq = pd.DataFrame(dataq)
    >>> check_context(dft, 'target', chrom_dict)
    Warning: added context to target is out of border for 3 instances
       start_1st chrom_1st  end_1st
    0          0      chr1       30

    >>> check_context(dfq, 'query', chrom_dict)
    Warning: added context to query is out of border for 4 instances
       start_2end  end_2end chrom_2end
    0           0        30       chr1


        """
    #print(df.info())
    #print(seq_tag)
    no_seq_out_boder = 0
    if seq_tag == 'target':
        start = 'start_1st'
        end = 'end_1st'
        chrom = 'chrom_1st'
    elif seq_tag == 'query':
        start = 'start_2end'
        end = 'end_2end'
        chrom = 'chrom_2end'

    df_copy = df.copy()
    df_filted, count = filter_false_chr(df_copy, chrom)
    no_seq_out_boder += len(df_filted[df_filted[start] < 0])
    no_seq_out_boder += len(df_filted[df_filted[end] > df_filted[chrom].apply(lambda x: chrom_dict[x])])

        # delet data
    df_filted = df_filted.loc[df_filted[start] >= 0]
    df_filted = df_filted.loc[df_filted[end] <= df_filted[chrom].apply(lambda x: chrom_dict[x])]

    if no_seq_out_boder > 0 :
        print(f'Warning: added context to {seq_tag} is out of border for {no_seq_out_boder} instances')
    return df_filted, count


################################################################################

def filter_false_chr(df, col_name):
    """
    If a cell of a column has False as a entry the row will be filtered out!

        Parameters
        ----------
        df: DataFrame
        col_name: chromosome column name


        Returns
        -------
        df_filtered
            dataframe without rows containing False rows with the col_name column
            df with changed positions

    >>> data = {'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'False', 'chr1', 'False']}
    >>> df = pd.DataFrame(data)
    >>> filter_false_chr(df, 'chrom_1st')
    (   start_1st chrom_1st
    0          0      chr1
    2         40      chr1, 2)
        """

    df_filterd = df[(df[col_name] != False)&(df[col_name] != 'False')]
    no_del_entys = len(df) - len(df_filterd)
    return df_filterd, no_del_entys


################################################################################

def get_context(seq_tag, df, out_dir, in_genome_fasta, context, chrom_len_file):
    """
    defining column with ID and empty columns to store the context sequences

        Parameters
        ----------
        seq_tag: dataframe
        df: dataframe containing position of the extraction
        out_dir: directory where to store bed and fa file
        in_genome_fasta: genome FASTA file
        context: amount of nt that should be added on both sites
        chrom_len_file: length of chromosomes file (Tabular)


        Returns
        -------
        df
            column update dataframe

        """
    out_bed = out_dir + seq_tag + '_out.bed'
    out_fa = out_dir + seq_tag + '_out.fa'
    no_del_entys = 0
    chrom_dict = read_table_into_dic(chrom_len_file)
    #print(chrom_dict)
    if seq_tag == 'target':
        df_bed = df[['chrom_1st', 'start_1st', 'end_1st', 'ID1', 'interaction_no', 'strand_1st']].copy()
        # print('test why False!!')
        df_bed['chrom_1st'] = df_bed['chrom_1st'].apply(lambda x: check_convert_chr_id(x))
        # df_bed.to_csv(out_bed, sep="\t", index=False)
        #df_bed_filted, count2 = filter_false_chr(df_bed, 'chrom_1st')
        #print(len(df_bed_filted))
        df_context =  add_context(df_bed, context, 'start_1st', 'end_1st')
        #print(df_context.tail())
        df_context, count = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_target'
        col_id = 'ID1'
        # df_context_filted, count = filter_false_chr(df_context, 'chrom_1st')
        no_del_entys += count
        #df_context_filted = df_context[df_context.chrom_1st != False]
        #no_del_entys += len(df_context) - len(df_context_filted)
    elif seq_tag == 'query':
        df_bed = df[['chrom_2end', 'start_2end', 'end_2end', 'ID2', 'interaction_no', 'strand_2end']].copy()
        df_bed['chrom_2end'] = df_bed['chrom_2end'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_2end', 'end_2end')
        df_context, count = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_query'
        col_id = 'ID2'
        # df_context_filted, count = filter_false_chr(df_context, 'chrom_2end')
        no_del_entys += count
        #df_context_filted = df_context[df_context.chrom_2end != False]
        #no_del_entys += len(df_context) - len(df_context_filted)
    else:
        print('error: please specify the parameter seq_tag with target or query')
    # delet all 'False' chromosomes of in the df
    print('lost %i instances because of the chromosome'%(no_del_entys))
    df_context.to_csv(out_bed, sep="\t", index=False, header=False)
    #df = df_context
    seqs_dic = bed_extract_sequences_from_fasta(out_bed, out_fa, in_genome_fasta)

    for seq_id in seqs_dic:
        #print(seq_id)
        #print(seqs_dic[seq_id])
        df.loc[df[col_id] == seq_id, [col_name]] = seqs_dic[seq_id]

    return df


################################################################################
#Functions negative data

def shuffle_sequence(seq, times, kind_of_shuffel):
    """
    shuffle on given sequence x times

        Parameters
        ----------
        seq: sequence
        times: amount of shuffling
        kind_of_shuffel: 1 -> Mononucleotide; 2 -> Dinucleotide


        Returns
        -------
        seq_list
            list of shuffled sequences

        """
    #seq_list= []

    call = "ushuffle -s " + str(seq) + " -n " \
                       + str(times) + " -k " + str(kind_of_shuffel)
    # print(call)
    p = subprocess.Popen(call, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True,)
    stdout, stderr = p.communicate()
    seq_list = stdout.decode('utf-8').strip().split('\n')

    return seq_list


################################################################################

def bp_suffeling(hybrid_seq, IntaRNA_prediction,times):
    """
    base pair shuffling of the given IntaRNA prediction

        Parameters
        ----------
        hybrid_seq: target sequence
        IntaRNA_prediction: query sequence


        Returns
        -------
        target
            list of shuffled target sequences
        shuffled_query_list
            list of shuffled query sequences

        """
    shuffled_target_list =[]
    shuffled_query_list = []

    tup_list = encode_hybrid_by_BPs(IntaRNA_prediction, hybrid_seq)
    #print(tup_list)
    # randomize the list with tuples where each tuple is a bp or bulge
    for i in range(times):
        shuffled_list = random.sample(tup_list, k=len(tup_list))
        traget, query = make_seq_from_list(shuffled_list)
        shuffled_target_list.append(traget)
        shuffled_query_list.append(query)

    return shuffled_target_list, shuffled_query_list

################################################################################

def load_occupied_data(input_occupied):
    """
    load occupied data

        Parameters
        ----------
        input_occupied: input file path

        Returns
        -------
        occupied_InteLab: Interlab object of occupied regions
        """
    overlap_handle = open(input_occupied,'rb')
    occupied_InteLab = pickle.load(overlap_handle)
    # print(overlap_avg_val)
    overlap_handle.close()
    return occupied_InteLab

################################################################################

def encode_hybrid_by_BPs(dot_bracket, seq):
    """
    encode_hybrid_by_BPs

        Parameters
        ----------
        dot_bracket: interaction given in with a dot bracket encoding
        seq: query

        Returns
        -------
        tup_list
            List of tuples storing the complete interaction of both
            [(target,query)]

        """
    # Test:    seq = 'ACCCACCCCCAA&AAGGAAGGGGGGA' hybrid = '.(((.(((((..&..))..)))))).'
    # result: [('A', '-'), ('-', 'A'), ('C', 'G'), ('C', 'G'), ('C', 'G'), ('A', '-'), ('C', 'G'), ('C', 'G'), ('C', 'G'), ('-', 'A'), ('-', 'A'), ('C', 'G'), ('C', 'G'), ('A', '-'), ('-', 'A'), ('A', '-'), ('-', 'A')]
    dot_bracket_list = list(dot_bracket)
    seq_list = list(seq)

    assert len(dot_bracket_list) == len(seq_list), 'RRI sequence and dot bracket string do not have the same length'

    idx_end = len(seq_list) - 1
    idx_start = 0
    tup_list = []
    for idx, start in enumerate(dot_bracket_list):
        end = dot_bracket_list[idx_end]
        start = dot_bracket_list[idx_start]
        if start == '&' and end == '&':
            break
        elif start == '(' and end == ')':
            tup_list.append((seq_list[idx_start],seq_list[idx_end]))
            idx_end -= 1
            idx_start += 1
        elif start == '.' and end == '.':
            tup_list.append((seq_list[idx_start],'-'))
            tup_list.append(('-',seq_list[idx_end]))
            idx_start += 1
            idx_end -= 1
        elif start == '.':
            tup_list.append((seq_list[idx_start],'-'))
            idx_start += 1
        elif end == '.':
            tup_list.append(('-',seq_list[idx_end]))
            idx_end -= 1
        else:
            print('hybrid encode error: unexpected case')
    return tup_list


################################################################################

def make_seq_from_list(shuffled_list):
    """
    make_seq_from_list

        Parameters
        ----------
        shuffled_list: shuffled list


        Returns
        -------
        seq1
            shuffled target sequences
        seq2
            shuffled query sequence

        """
    seq1 = ''
    seq2 = ''
    for tup in shuffled_list:
        character_seq1 = tup[0]
        character_seq2 = tup[1]
        if character_seq1 != '-' and character_seq2 != '-':
            seq1 = seq1 + character_seq1
            seq2 = seq2 + character_seq2
        elif character_seq1 == '-':
            seq2 = seq2 + character_seq2
        elif character_seq2 == '-':
            seq1 = seq1 + character_seq1
        else:
            print('hybrid encode error: something went wrong with the encoding')

    return seq1, seq2


################################################################################

def mearge_overlaps(inter_obj, info):
    """
    merge positions in a interlab library

        Parameters
        ----------
        inter_obj : inter objects
        info: information what inter object


        Returns
        -------
        inter_obj_new
            interlap objects with the merged positions

        """
    inter_obj_new = defaultdict(InterLap)

    for key in inter_obj:
        #print(key)
        #print(list(inter_obj[key]))
        inter_list_temp = [(i[0],i[1]) for i in list(inter_obj[key])]
        #print(inter_list_temp)
        inter = join_pos(inter_list_temp)
        #print(inter)
        inter_list = [(i[0],i[1], info) for i in list(inter)]
        #print(inter_list)
        inter_obj_new[key].add(inter_list)
        #for i in inter_rep_one[key]:
            #print(i)
        #print('test interval')
    return inter_obj_new


################################################################################

def join_pos(pos_list):
    """
    join positions will join start end end positions which are overlapping

        Parameters
        ----------
        pos_list : list of tuples containing (start, end) position

        Returns
        -------
        inter_obj_new
            interlap objects with the merged positions

    >>> join_pos([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]
    >>> join_pos([(2, 6), (4, 10)])
    [(2, 10)]
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


################################################################################

def read_table_into_dic(file):
    """
    Read in Table separated by \t and puts first line as key second as value
        Parameters
        ----------
        file: file location at the table file

        Returns
        -------
        chrom_ends_dic:
            chrom -> length

    """
    chrom_ends_dic = {}

    # Open FASTA either as .gz or as text file.
    if re.search(".+\.gz$", file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file, "r")

    for line in f:
        line_list = line.split("\t")
        chrom_ends_dic[line_list[0]] = int(line_list[1].rstrip())
    f.close()

    return chrom_ends_dic


################################################################################

def read_pos_neg_data(in_positive_data_filepath, in_negative_data_filepath):
    pos_df = pd.read_csv(in_positive_data_filepath, sep=',')
    neg_df = pd.read_csv(in_negative_data_filepath, sep=',')
    #Inject labels
    pos_df['label'] = 1
    neg_df['label'] = 0
    #Dataset initial characterization
    reporting=0
    if(reporting):
        pos_report=pandas_profiling.ProfileReport(pos_df,title="Positive data Report")
        neg_report=pandas_profiling.ProfileReport(neg_df,title="Negative data Report")
        pos_report.to_file(output_path + "/positive_report.html")
        neg_report.to_file(output_path + "/negative_report.html")
    #print(pos_df.dtypes)
    #print(neg_df.dtypes)
    #print(pd.get_dummies(pos_df))
    #print(pd.get_dummies(neg_df))
    #Concat datasets
    ia_df = pd.concat([pos_df,neg_df])

    y = ia_df.label
    X = ia_df.drop(columns="label")

    return X, y


################################################################################
#Functions for model training

def train_model(in_positive_data_filepath,in_negative_data_filepath,output_path):

    X, y = read_pos_neg_data(in_positive_data_filepath, in_negative_data_filepath)

    #y = ia_df.label
    #X = ia_df.drop(columns="label")
    for m in [DummyClassifier, LogisticRegression, DecisionTreeClassifier, KNeighborsClassifier,GaussianNB, SVC, RandomForestClassifier, xgboost.XGBClassifier]:
        cls=m()
        kfold = model_selection.KFold(n_splits=10, random_state=42, shuffle=True)
        s = model_selection.cross_val_score(cls, X,y, scoring="roc_auc", cv=kfold)
        print(
            f"{m.__name__:22}\t AUC:\t"
            f"{s.mean():.3f}\t STD:\t {s.std():.2f}"
            )
    #Create training and test dataset
    X_training, X_test, y_training, y_test = model_selection.train_test_split(X, y, test_size=0.3, random_state=42)
    ##comparison dummy model
    cm = DummyClassifier()
    cm.fit(X_training, y_training)
    dummy_comparison_score = cm.score(X_test, y_test)
    print("Dummy score: %f" %(dummy_comparison_score))
    #random_forest
    random_forest = RandomForestClassifier(n_estimators=100, random_state=42)
    random_forest.fit(X_training, y_training)
    random_forest_comparison_score = random_forest.score(X_test, y_test)
    print("RF score: %f" %random_forest_comparison_score)
    rf_path = output_path + "/rf.obj"
    rf_handle = open(rf_path,"wb")
    pickle.dump(random_forest,rf_handle)
    rf_handle.close()

    xgb = xgboost.XGBClassifier(n_estimators=100, random_state=42)
    xgb.fit(X_training, y_training)
    xgb_comparison_score = xgb.score(X_test, y_test)
    print("RF score: %f" %xgb_comparison_score)
    xgb_path = output_path + "/xgb.obj"
    xgb_handle = open(xgb_path,"wb")
    pickle.dump(xgb,xgb_handle)
    xgb_handle.close()
    return ""


################################################################################

def base_model(in_positive_data_filepath,in_negative_data_filepath,output_path,name):

    X, y = read_pos_neg_data(in_positive_data_filepath, in_negative_data_filepath)

    #Create training and test dataset
    X_training, X_test, y_training, y_test = model_selection.train_test_split(X, y, test_size=0.3, random_state=42)
    ##comparison dummy model
    cm = DummyClassifier()
    cm.fit(X_training, y_training)
    # prediction: “prior”: always predicts the class that maximizes the class prior (like “most_frequent”) and predict_proba returns the class prior.
    params = cm.get_params(deep=True)
    dummy_comparison_score = cm.score(X_test, y_test)
    print("Dummy score: %f" %(dummy_comparison_score))
    # save
    dc_path = output_path + "/" + name + "_dc.obj"
    dc_handle = open(dc_path,"wb")
    pickle.dump(cm,dc_handle)
    dc_handle.close()

    # C-Support Vector Classification
    svc = SVC()
    svc.fit(X_training, y_training)
    # mean accuracy
    svc_comparison_score = svc.score(X_test, y_test)
    print("Support Vector Classification score: %f" %svc_comparison_score)
    svc_path = output_path + "/" + name + "_svc.obj"
    svc_handle = open(svc_path,"wb")
    pickle.dump(svc,svc_handle)
    svc_handle.close()

    #print('trainings data X:\n', X_training.info())
    #print('trainings data y:\n', y_training.dtype)

    # linear regerssion -> LIBLINEAR – A Library for Large Linear Classification
    clf = LogisticRegression(solver='liblinear')
    clf.fit(X_training, y_training)

    logistic_regression_score = clf.score(X_test, y_test)
    print("LogisticRegression score: %f" %(logistic_regression_score))
    clf_path = output_path + "/" + name + "_clf.obj"
    clf_handle = open(clf_path,"wb")
    pickle.dump(clf,clf_handle)
    clf_handle.close()

    return ""


################################################################################


def get_ids(pos_data):
    """
    Extract target und query IDs from df
        ----------
        pos_data: file contining positive IntaRNA calls including ID's

        Returns
        -------
        df_ID
            datafram contianing IDs

    """
    df_pos = pd.read_table(pos_data, sep=',')
    df_ID = df_pos[['target_ID','query_ID']]

    return df_ID



def classify(df_eval,in_model_filepath, output_path, df_ID='off',
            true_lable=False, y='off', ID_Label = False):
    """
    Classification of the given RRIs features using the input model
        Parameters
        ----------
        df_eval: dataframe with interaction features of to evaluate instances
        in_model_filepath: path to the model
        output_path: location where to save output
        df_ID: all ids adding to the final table off for cross evaluation
        true_lable: Boolien to set lable or not (needed for evaluation)
        y: list containing labes to the df_eval instaces
        ID_Label: check if ID should be added


        Returns
        -------
        df_result
            data frame haveing all columns of df_eval and prediction column
            this df is also saved as csv in output

    """
    #X = pd.read_csv(in_data_filepath, sep=',')
    #model_handle = open(in_model_filepath,'rb')
    #model = pickle.load(model_handle)
    #model_handle.close()
    #params = loadfile(in_model_filepath)['params']
    #print(params)
    model = loadfile(in_model_filepath)['estimator']
    #model = loadfile(in_model_filepath)
    y_pred=model.predict(df_eval)
    #print('model predictions')
    xtra = pd.DataFrame({'predicted_label': y_pred})
    df_result = pd.concat([xtra, df_eval], axis=1)
    # df_eval['prediction'] = y_pred


    if true_lable:
        df_result['true_label'] = y.tolist()
        # calculates prob for each label. We only returen prob for one label?
    proba = model.predict_proba(df_eval)[:,1]
    instance_score = pd.DataFrame({'instance_score': proba})
    #print(instance_score)
    #print(df_result)
    if ID_Label:
        df_result_temp = pd.concat([df_ID, instance_score, df_result], axis=1)
        df_result = df_result_temp.loc[df_result_temp.groupby(['target_ID', 'query_ID'])['instance_score'].idxmax()]
        # implemnt a majorety vote
        # df_data.groupby(['target_ID', 'query_ID']).agg({'predicted_label': lambda x: x.value_counts().index[0]}).reset_index()
    else:
        df_result = pd.concat([instance_score, df_result], axis=1)

    # reformat data structure

    df_result.to_csv(output_path,index=False)
    return df_result


def extend_eval_with_empty_entries(df_eval, pos_neg_out_path, output_path, substract_numbers=4):
    """
    Extend the final output table with all instances that had no IntaRNA 
    prediction and can therefore not be predicted.

        Parameters
        ----------
        df_eval: dataframe with interaction features of to evaluate instances
        pos_neg_out_path: path to positive neg data where the json file is stored
        output_path: location where to save the final evaluation table
        substract_numbers: number of additional colums that are not feature  
                           (instance_score, prediction, ...)
        

    """
    num_feat_cols = len(df_eval.columns) - substract_numbers  

    na_rows_feat = [np.nan] * num_feat_cols
    #print(len(df_eval))
    header = df_eval.columns

    # read jason file 
    with open(f'{pos_neg_out_path}/missed_list.json', 'r') as file:
        id_list = json.load(file)
    extend_list = []
    print(f'Length of missed IDs: {len(id_list)}')
    for id in id_list:
        # print(type(id))
        new_row = []
        id_temp = id.split('_')
        # print(id_temp)
        new_row = [id_temp[0], id_temp[1], 0, 0] + na_rows_feat
        print(f'new row: {new_row}')
        extend_list.append(new_row)
        #print(len(df_eval))
        #df_eval.loc[(len(df_eval)+1)] = new_row
        # df_eval = df_eval.append(new_row, ignore_index=True)
    print(df_eval)
    df_extention = pd.DataFrame(extend_list, columns=df_eval.columns)
    print(df_extention)   
    df_concat = pd.concat([df_eval, df_extention], ignore_index=True)

    df_concat.to_csv(output_path,index=False)



################################################################################

def classify2(X,in_model_filepath,score_flag):
    model_handle = open(in_model_filepath,'rb')
    model = pickle.load(model_handle)
    model_handle.close()
    y_pred = model.predict(X)
    if score_flag == 'proba':
        y_score = model.predict_proba(X)
    elif score_flag == 'decision':
        y_score = model.decision_function(X)
    #print('model predictions')
    #print(y_pred)
    return model, y_pred, y_score


################################################################################

def param_optimize(in_positive_data_filepath,in_negative_data_filepath,output_path):
    X, y = read_pos_neg_data(in_positive_data_filepath, in_negative_data_filepath)

    #y = ia_df.label
    #X = ia_df.drop(columns="label")
    X_training, X_test, y_training, y_test = model_selection.train_test_split(X, y, test_size=0.3, random_state=42)

    # computing base model performance:
    base_model = RandomForestClassifier(n_estimators=100, random_state=42)
    base_model.fit(X_training, y_training)
    base_accuracy = evaluate(base_model, X_test, y_test)


    #random_forest
    random_forest = RandomForestClassifier()
    # dict of hyperparameters to optimize
    param_grid = {'bootstrap': [True],
        'max_depth': [6, 10],
        'max_features': ['auto', 'sqrt'],
        'min_samples_leaf': [3, 5],
        'min_samples_split': [4, 6],
        'n_estimators': [100, 350]
        }

    forest_grid_search = GridSearchCV(random_forest, param_grid, cv=5,
                                      scoring="roc_auc",
                                      return_train_score=True,
                                      verbose=True,
                                      n_jobs=-1)

    forest_grid_search.fit(X_training, y_training)

    best_param = forest_grid_search.best_params_
    print("RF best params: ")
    print(best_param)

    best_grid = forest_grid_search.best_estimator_
    grid_accuracy = evaluate(best_grid, X_test, y_test)

    print('RF base accuracy: %f' % base_accuracy)
    print('RF base accuracy: %f' % grid_accuracy)

    print('Improvement of {:0.2f}%.'.format( 100 * (grid_accuracy - base_accuracy) / base_accuracy))



    #random_forest_comparison_score = random_forest.score(X_test, y_test)
    #print("RF score: %f" %random_forest_comparison_score)
    #rf_path = output_path + "/rf.obj"
    #rf_handle = open(rf_path,"wb")
    #pickle.dump(random_forest,rf_handle)
    #rf_handle.close()

    return ""


################################################################################

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))

    return accuracy


################################################################################

def filter_features(X,featurefile,use_structure):
    if use_structure == 'on':
        # feature selection is called within the graph feature generation
        ft = np.load(featurefile)['d']
    elif use_structure == 'off':
        # file is generated by calling biofilms feature selection
        ft = np.load(featurefile)['arr_0']
    #print(X.info())
    header = X.columns
    header_list = header.tolist()
    features_filtered=list(compress(header_list, ft))
    #print(features_filtered)
    X_filtered = X[features_filtered]
    #print('dfInfo after:')
    #print(X_filtered.info())

    return X_filtered



################################################################################



################################################################################

def get_filted_features(featurefile,csr_mat_h):
    """
    produces a Boolean header based on a already filtered real header so the
    evaluation data can filter for the same columns as was used in training
        Parameters
        ----------
        featurefile: feature file used for model training


        Returns
        -------
        feat
            Boolean list to filter evaluation feature data

    """
    ft = np.load(featurefile)['d']

    #ft_new = np.append('index', ft)
    #ft_new = np.insert(ft, 0, 'index')
    header = ['E', 'E_hybrid', 'ED1', 'ED2',
                'len_interaction_target', 'len_interaction_query', 'no_bps',
                'max_inter_len', 'inter_len_normby_bp',
                'bp_normby_inter_len', 'GC_content', 'GC_skew', 'AT_skew',
                'max_ED', 'sum_ED', 'mfe_normby_GC', 'max_ED_normby_GC',
                'E_hybrid_normby_GC', 'mfe_normby_len', 'max_ED_normby_len',
                'E_hybrid_normby_len', 'mfe_normby_GC_len',
                'max_ED_normby_GC_len', 'E_hybrid_normby_GC_len',
                'complex_target_site', 'complex_query_site']
    dataset_size = len(csr_mat_h)
    feat = [False] * dataset_size

    for el in ft:
        # get pos if handcrafted featurs
        if el in header:
            pos = header.index(el)
            feat[pos] = True
        else:
            #el_int = el.astype('int')
            #print(int(el))
            feat[(int(el)+len(header))] = True

    return feat



def write_json_index(list_dfs, file_json):
    """
    takes a list of dataframes concatintes them and builds a index -> ID
    dictionary which is saved into a json file. The df's should have the same
    clolums with one of it beeing the ID
    (e.g. chr15;+;82519671;82519691_chrX;+;109054541;109.)
        Parameters
        ----------
        list_dfs: list of dataframes e.g. [df_pos,df_neg]
        file_json: file pathe to the output json file


    """
    data_df = pd.concat(list_dfs, ignore_index=True)
    #print(data_df)
    #print(len(data_df))
    data_ID_dict = data_df['ID'].to_dict()

    # print(len(data_ID_dict))
    write_json_file(data_ID_dict, file_json)
    #with open(file_json, 'w') as file:
        #json.dump(data_ID_dict, file)


def write_json_file(content, file_json):
    """
    takes a python object and dumps it into a json file
        Parameters
        ----------
        content: python object
        file_json: file pathe to the output json file


    """
    with open(file_json, 'w') as file:
        json.dump(content, file)



################################################################################


def call_vectorize(g):
    return eg.vectorize(g,nbits=18)

def convert(X, y, outname, graphfeatures, mode, feat_file='non', no_jobs=1):

    #call_script(f'export PYTHONHASHSEED=31337')

    # makes list of subseqDP and hybridDP tupels
    hybrid_seq_list = [a for a in zip(X['subseqDP'],X['hybridDP'])]
    X = X.drop(columns="subseqDP")
    X = X.drop(columns="hybridDP")

    #print(X.index.tolist())
    # X['index'] = X.index.tolist()
    h_handcrafted_feat = X.columns.tolist()
    y_np = np.array(y.tolist())
    #print(h_handcrafted_feat)

    if graphfeatures:
        # convert df into a csr matrix
        # print(f'number of jobs: {no_jobs}')
        X_from_df = csr_matrix(X.to_numpy().astype(np.float64))
        graphs = tools.xmap(mkgr, hybrid_seq_list, no_jobs)

        graphs_csr = csr_matrix(vstack(tools.xmap(call_vectorize,[[g] for g in graphs],no_jobs)))
        no_graphs = graphs_csr.get_shape()[1]
        #print('computed call_verctorize!!')
        X_csr= csr_matrix(hstack((X_from_df,graphs_csr)))

        # add eden features to header
        h_eden_feat = [str(int) for int in np.arange(0, no_graphs, 1).tolist()]
        header_full = h_handcrafted_feat + h_eden_feat

        # call feature selection
        if mode == 'train':
            feat,_ = featsel.forest(X_csr,y_np,X_csr,y_np,None)
        elif mode == 'eval':
            #feat = np.load(feat_file)['arr_0']
            feat = get_filted_features(feat_file,header_full)
        X_csr = X_csr[:,feat]
        header = [d for d, s in zip(header_full, feat) if s]

        #print(header)

        X_np= X_csr.todense()

    else:
        X_np = X.to_numpy()
        header = h_handcrafted_feat

    # saves object as 'np.savez_compressed'
    index = X.index.tolist()
    list = [X_np, y_np, header, index]
    tools.ndumpfile(list , outname)
    #print(f'total amount of features:\n {len(X.columns.tolist())}')
    #print()
    if mode == 'eval':
        #print(header)
        X = pd.DataFrame(X_np, columns=header)
        #print(header)
        #X['true_lable'] = y
        #print(len(y.tolist()))
        #print(X.info())
        return X,y



def mkgraph(sequence, structure):
    graph = nx.Graph()
    lifo = defaultdict(list)
    cut = structure.index("&")

    structure = structure.replace("&","")
    sequence = sequence.replace("&","")
    for i,(s,n) in enumerate(zip(structure, sequence)):
        graph.add_node(i, label=n)
        if i > 0 and  i != cut:
            graph.add_edge(i, i-1, label='-')

        # ADD PAIRED BASES
        if s in ['(','[','<']:
            lifo[partner[s]].append(i)
        if s in [')',']','>']:
            j = lifo[s].pop()
            graph.add_edge(i, j, label='=')
    return graph
    #return eg.vectorize([graph], discrete = False) # keep here in case i want nested edges ...


def mkgr(x):
    return mkgraph(*x)

def download_file(url, outdir):
    filename = wget.download(url, out=outdir)


def download_genome(out_path, genome, chrom_len_file):
    genome_dir = f'{out_path}/genome/'
    if not os.path.exists(genome_dir):
        os.mkdir(genome_dir)
    if genome == 'human':
        genome_file = f'{genome_dir}hg38.fa'
        chrom_len_file = f'{genome_dir}hg38.chrom.sizes'
        if not os.path.isfile(genome_file):
            genome_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz'
            download_file(genome_url, genome_dir)
            os.system("gunzip " + genome_file + ".gz")
        if not os.path.isfile(chrom_len_file):
            chom_len_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
            download_file(chom_len_url, genome_dir)
    elif genome == 'mouse':
        genome_file = f'{genome_dir}mm10.fa'
        chrom_len_file = f'{genome_dir}mm10.chrom.sizes'
        if not os.path.isfile(genome_file):
            genome_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz'
            download_file(genome_url, genome_dir)
            os.system("gunzip " + genome_file + ".gz")
        if not os.path.isfile(chrom_len_file):
            chom_len_url = 'https://hgdownload.soe.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes'
            download_file(chom_len_url, genome_dir)
    else:
        if not os.path.isfile(genome):
            print(f'Error: please provide genome FASTA file path')
        else:
            genome_file = genome
        if not os.path.isfile(chrom_len_file):
            print(f'Error: please provide a chromosome length file path')
        else:
            chrom_len_file = chrom_len_file
    return (genome_file, chrom_len_file)


#################################################################################
def get_files_concat(file_path,out_file_prefix):
    """
    From a list of files concatinate all pos/neg files given in the list, saves
    the the cocatinated pos/neg files
        Parameters
        ----------
        file_path: list containing filepathes the to be concatinated dataq
        out_file_prefix: prefix output file path

        Returns
        -------
        pos_file_out
            file location of the postive data file
        neg_file_out
            file location of the negative data file

    """
    pos_list = []
    neg_list = []

    for file_path in file_path:

            #feature_path = (f'{input_path_RRIs}/{data}/feature_files/'
            #               f'feature_filtered_{data}_context_{str(context)}_pos_occ_')

        file_neg = (f'{file_path}neg.csv')
        file_pos = (f'{file_path}pos.csv')


            #X_list.append(X_sub)
            #y_list.append(y_sub)

        pos_list.append(pd.read_csv(file_pos, sep=','))
        neg_list.append(pd.read_csv(file_neg, sep=','))

        #X = pd.concat(X_list)
        #y = pd.concat(y_list)

    df_pos = pd.concat(pos_list)
    df_neg = pd.concat(neg_list)
    pos_file_out = (f'{out_file_prefix}_pos_occ_pos.csv')
    neg_file_out = (f'{out_file_prefix}_pos_occ_neg.csv')

    df_pos.to_csv(pos_file_out,sep=",", index=False)
    df_neg.to_csv(neg_file_out,sep=",", index=False)

    return pos_file_out, neg_file_out


def perfome_cv(folds, opt_call_cv, out_path_model, midel_name):
    """
    Cross validation for the model perfomance
        Parameters
        ----------
        folds: number of folds that should be perfomed
        opt_call_cv: prefix optimization calls
        out_path_model: pefix of the output model file name
        midel_name: exp name and context string to build file name

        Returns
        -------
        f1
            F1-score

    """
    for fold in range(folds):
        cv_call = (f'{opt_call_cv} --folds {folds} --foldselected {fold}'
                   f' --out {out_path_model}{midel_name}_fold{fold} '
                   f'--autosk_debugfile {out_path_model}{midel_name}_'
                   f'fold{fold}_autosklearn.log '
                   f'--autosk_debug True')
        print(f'\n run optimize call for fold {fold}')
        # print(cv_call)
        out = call_script(cv_call, reprot_stdout=True, asset_err=False)
    f1_df_list = []
    for fold in range(folds):
        filename = f'{out_path_model}{midel_name}_fold{fold}.csv'
        df_temp = pd.read_csv(filename)
        f1_df_list.append(df_temp)

    df_f1 = pd.concat(f1_df_list, axis=0)
    # print(len(df_f1))
    f1 = f1_score(df_f1['true_label'].tolist(), df_f1['predicted_label'].tolist())
    return f1
