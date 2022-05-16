#!/usr/bin/env python

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



"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
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
        in_file : tabular file with output of chira RRI results

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
                'ineraction_site_1st', 'ineraction_site_2end',
                'IntaRNA_prediction', 'energy',
                'seq_1st_ineraction_site', 'seq_2end_ineraction_site',
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

def filter_score(df, score_th):
    """
    Filter dataframe for instances with a score of 1

        Parameters
        ----------
        df : df including the containing all RRIs (Interactions)
        score_th: threshold for the expactation maximization score of Chira


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
    Filters out all rows with a emty column entry in a dataframe

        Parameters
        ----------
        df : df including the containing all RRIs
        col_name: name of to be filterd colum


        Returns
        -------
        df_filtered
            dataframes without empty column entrys


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
    Starts a subprosses to call a script and checks for errors.


        Parameters
        ----------
        call : cmd comand

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
    Building for each replicat a inter object

        Parameters
        ----------
        s1: start of one sequence of the first replicat
        e1: end of one sequence of the first replicat
        s2: start of one sequence of the current replicat
        e2: end of one sequence of the current replicat

        Returns
        -------
        compinde_overlap
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
            print('error: somthing is not overlaping hier')
    elif s2 < s1:
        s_overlap = s1
        if e1 <= e2:
            e_overlap = e1
        elif e2 < e1:
            e_overlap = e2
        else:
            print('error: somthing is not overlaping hier')
    overlap_len = e_overlap - s_overlap +1
    seq1_len = e1 - s1 + 1
    seq2_len = e2 - s2 + 1
    overlap_seq1 = overlap_len/seq1_len
    # print(overlap_seq1)
    overlap_seq2 = overlap_len/seq2_len
    # print(overlap_seq2)
    # compinde_overlap = (overlap_seq1 + overlap_seq2)/2
    # select overlap of shorter sequence:
    compinde_overlap = max([overlap_seq1, overlap_seq2])
    # print(compinde_overlap)
    if len_flag:
        return overlap_len
    else:
        return compinde_overlap


################################################################################

def get_chrom_list_no_numbers(df_interactions, chrom):
    """
    Generates a unique list of chromosmes or conticts for both interaction
    partners

        Parameters
        ----------
        df_interactions : df including the filtered RRIs
        chrom : string indicating from wich seq the chromosome is


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


################################################################################

def get_list_chrom(df):
    """
    Generates a unique list of chromosmes or conticts for both interaction
    partners

        Parameters
        ----------
        df : df including the filtered RRIs


        Returns
        -------
        sort_list_chrom
            sorted list of unique chromosmes or contics which are not a number
            and present in the input data frame

        """
    chrom1_list = get_chrom_list_no_numbers(df, 'chrom_seq_1st_site')
    chrom2_list = get_chrom_list_no_numbers(df, 'chrom_seq_2end_site')
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
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

        Parameters
        ----------
        chr_id: chromosme id string


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


################################################################################

def add_context(df_bed, context, start, end):
    """
    edding the changing the start and end postion of the sequences
    to add context to both sites of the sequences in the dataframe

        Parameters
        ----------
        df_bed: dataframe containing start and end positon
        context: amount of nucleotied
        start: column name of start positons
        end: column name of end positons


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


################################################################################

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
        seqs_dic = read_fasta_into_dic(out_fa, skip_n_seqs=True)
        # Output RNA sequences.
        #fasta_output_dic(seqs_dic, out_fa, split=True)
    return seqs_dic


################################################################################

def check_context(df, seq_tag, chrom_dict):
    """
    check that the extended context is not to short or long!

        Parameters
        ----------
        df: df with added context
        seq_tag: target or query
        chrom_dict: dictionary storing chrmosome length {name->len,...}


        Returns
        -------
        df
            df with changed postions

    >>> chrom_dict = {'chr1':60,
    ...         'chr2':80}
    >>> datat = {'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'chr1', 'chr1', 'chr1'],
    ...         'end_1st':[30, 60, 70, 70]}
    >>> dataq = {'start_2end':[0, -40, 40, -2],
    ...         'end_2end':[30, 60, 70, 70],
    ...         'chrom_2end':['chr1', 'chr1', 'chr1', 'chr1']}
    >>> dft = pd.DataFrame(datat)
    >>> dfq = pd.DataFrame(dataq)
    >>> check_context(dft, 'target', chrom_dict)
    Warning: added context to target is out of bourder for 4 instances
       start_1st chrom_1st  end_1st
    0          0      chr1       30

    >>> check_context(dfq, 'query', chrom_dict)
    Warning: added context to query is out of bourder for 4 instances
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

    #print(df)
    no_seq_out_boder += len(df[df[start] < 0])
    no_seq_out_boder += len(df[df[end] > df[chrom].apply(lambda x: chrom_dict[x])])

        # delet data
    df = df.loc[df[start] >= 0]
    df = df.loc[df[end] <= df[chrom].apply(lambda x: chrom_dict[x])]

    print('Warning: added context to %s is out of bourder for %i instances'%(seq_tag,no_seq_out_boder))
    return df


################################################################################

def filter_false_chr(df, col_name):
    """
    If a cell of a column has False as a entry the row will be filterd out!

        Parameters
        ----------
        df: DataFrame
        col_name: chromosome column name


        Returns
        -------
        df_filtered
            dataframe without rows contining False rows with the col_name column
            df with changed postions

    >>> data = {'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'False', 'chr1', 'False']}
    >>> df = pd.DataFrame(data)
    >>> filter_false_chr(df, 'chrom_1st')
    (   start_1st chrom_1st
    0          0      chr1
    2         40      chr1, 2)
        """
    df_filterd = df[df[col_name] != 'False']
    no_del_entys = len(df) - len(df_filterd)
    return df_filterd, no_del_entys


################################################################################

def get_context(seq_tag, df, out_dir, in_2bit_file, context, chrom_len_file):
    """
    defining column with ID and empty colums to store the context sequences

        Parameters
        ----------
        seq_tag: dataframe
        df: dataframe contining position of the extraction
        out_dir: directory where to store bed and fa file
        in_2bit_file: genome 2bit file
        context: amout of nt that should be added on both sites


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
        #print(df_bed.tail())
        df_bed['chrom_1st'] = df_bed['chrom_1st'].apply(lambda x: check_convert_chr_id(x))
        #print(df_bed.tail())
        df_context =  add_context(df_bed, context, 'start_1st', 'end_1st')
        #print(df_context.tail())
        df_context = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_target'
        col_id = 'ID1'
        df_context_filted, count = filter_false_chr(df_context, 'chrom_1st')
        no_del_entys += count
        #df_context_filted = df_context[df_context.chrom_1st != False]
        #no_del_entys += len(df_context) - len(df_context_filted)
    elif seq_tag == 'query':
        df_bed = df[['chrom_2end', 'start_2end', 'end_2end', 'ID2', 'interaction_no', 'strand_2end']].copy()
        df_bed['chrom_2end'] = df_bed['chrom_2end'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_2end', 'end_2end')
        df_context = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_query'
        col_id = 'ID2'
        df_context_filted, count = filter_false_chr(df_context, 'chrom_2end')
        no_del_entys += count
        #df_context_filted = df_context[df_context.chrom_2end != False]
        #no_del_entys += len(df_context) - len(df_context_filted)
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
    basepair shufelling of the given IntaRNA prediction

        Parameters
        ----------
        hybrid_seq: tartet sequence
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
        suffled_list = random.sample(tup_list, k=len(tup_list))
        traget, query = make_seq_from_list(suffled_list)
        shuffled_target_list.append(traget)
        shuffled_query_list.append(query)

    return shuffled_target_list, shuffled_query_list

################################################################################

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

################################################################################

def encode_hybrid_by_BPs(dot_bracked, seq):
    """
    encode_hybrid_by_BPs

        Parameters
        ----------
        dot_bracked:
        seq: query

        Returns
        -------
        tup_list

        """
    # Test:    seq = 'ACCCACCCCCAA&AAGGAAGGGGGGA' hybrid = '.(((.(((((..&..))..)))))).'
    # result: [('A', '-'), ('-', 'A'), ('C', 'G'), ('C', 'G'), ('C', 'G'), ('A', '-'), ('C', 'G'), ('C', 'G'), ('C', 'G'), ('-', 'A'), ('-', 'A'), ('C', 'G'), ('C', 'G'), ('A', '-'), ('-', 'A'), ('A', '-'), ('-', 'A')]
    dot_bracked_list = list(dot_bracked)
    seq_list = list(seq)

    assert len(dot_bracked_list) == len(seq_list), 'RRI sequence and dotbracked string do not have the same lenght'

    idx_end = len(seq_list) - 1
    idx_start = 0
    tup_list = []
    for idx, start in enumerate(dot_bracked_list):
        end = dot_bracked_list[idx_end]
        start = dot_bracked_list[idx_start]
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
            print('hybrid encode error: unexpacted case')
    return tup_list


################################################################################

def make_seq_from_list(suffled_list):
    """
    make_seq_from_list

        Parameters
        ----------
        suffled_list:


        Returns
        -------
        seq1
            shuffled target sequences
        seq_2
            shuffled query sequence

        """
    seq1 = ''
    seq2 = ''
    for tup in suffled_list:
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
            print('hybrid encode error: soemthing went wrong with the encoding')

    return seq1, seq2


################################################################################

def mearge_overlaps(inter_obj, info):
    """
    mearg postions in a interlab library

        Parameters
        ----------
        inter_obj : inter objects
        info: information what inter object


        Returns
        -------
        inter_obj_new
            inerlap objects with the mearged positons

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
        file: file location ot the table file

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

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        dna=False,
                        report=1,
                        all_uc=False,
                        skip_data_id="set",
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    FASTA file can be plain text or gzipped (watch out for .gz ending).
        Parameters
        ----------
        fasta_file: file location of the to be read fasta file


        Returns
        -------
        seqs_dic
            dictonary with seq id as key and sequence as value

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""

    # Open FASTA either as .gz or as text file.
    if re.search(".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else:
        f = open(fasta_file, "r")
    for line in f:
        if re.search(">.+", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)
            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
        elif re.search("[ACGTUN]+", line, re.I):
            m = re.search("([ACGTUN]+)", line, re.I)
            seq = m.group(1)
            if seq_id in seqs_dic:
                if dna:
                    # Convert to DNA, concatenate sequence.
                    seq = seq.replace("U","T").replace("u","t")
                else:
                    # Convert to RNA, concatenate sequence.
                    seq = seq.replace("T","U").replace("t","u")
                if all_uc:
                    seq = seq.upper()
                seqs_dic[seq_id] += seq
    f.close()

    # Check if sequences read in.
    assert seqs_dic, "no sequences read in (input FASTA file \"%s\" empty or mal-formatted?)" %(fasta_file)
    # If sequences with N nucleotides should be skipped.
    c_skipped_n_ids = 0
    if skip_n_seqs:
        del_ids = []
        for seq_id in seqs_dic:
            seq = seqs_dic[seq_id]
            if re.search("N", seq, re.I):
                if report == 1:
                    print ("WARNING: sequence with seq_id \"%s\" in file \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id, fasta_file))
                c_skipped_n_ids += 1
                del_ids.append(seq_id)
        for seq_id in del_ids:
            del seqs_dic[seq_id]
        assert seqs_dic, "no sequences remaining after deleting N containing sequences (input FASTA file \"%s\")" %(fasta_file)
        if c_skipped_n_ids:
            if report == 2:
                print("# of N-containing %s regions discarded:  %i" %(skip_data_id, c_skipped_n_ids))
    return seqs_dic


################################################################################

def read_pos_neg_data(in_positive_data_filepath, in_negative_data_filepath):
    pos_df = pd.read_csv(in_positive_data_filepath, sep=',')
    neg_df = pd.read_csv(in_negative_data_filepath, sep=',')
    #Inject labels
    pos_df['label'] = 1
    neg_df['label'] = 0
    #Dataset initial characterisation
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

def classify(df_eval,in_model_filepath, output_path,true_lable=False, y='off'):
    """
    Classification of a given given elements using the input model
        Parameters
        ----------
        df_eval: dataframe with interaction features of to evalute instances
        in_model_filepath: path to the model
        output_path: location whter to save output

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
    #model = loadfile(in_model_filepath)['estimator']
    model = loadfile(in_model_filepath)
    y_pred=model.predict(df_eval)
    #print('model predictions')
    xtra = pd.DataFrame({'predicted_label': y_pred})
    df_result = pd.concat([df_eval, xtra], axis=1)
    # df_eval['prediction'] = y_pred
    #print(y_pred)
    if true_lable:
        df_result['true_label'] = y.tolist()
    df_result.to_csv(output_path,index=False)
    return df_result


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

    # computing base modle perfomance:
    base_model = RandomForestClassifier(n_estimators=100, random_state=42)
    base_model.fit(X_training, y_training)
    base_accuracy = evaluate(base_model, X_test, y_test)


    #random_forest
    random_forest = RandomForestClassifier()
    # dict of hyperparmeters to optimize
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
    features_filterd=list(compress(header_list, ft))
    #print(features_filterd)
    X_filterd = X[features_filterd]
    #print('dfInfo after:')
    #print(X_filterd.info())

    return X_filterd



################################################################################



################################################################################

def get_filted_features(featurefile,csr_mat_h):
    """
    produces a boolian header based on a alrady filtered real header so the
    evaluation data can filter for the same colums as was used in training
        Parameters
        ----------
        featurefile: feature file used for model training


        Returns
        -------
        feat
            boolian list to filter evaluation feature data

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
        print(f'number of jobs: {no_jobs}')
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
    list = [X_np, y_np, header]
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
        genome_file = f'{genome_dir}hg38.2bit'
        chrom_len_file = f'{genome_dir}hg38.chrom.sizes'
        if not os.path.isfile(genome_file):
            genome_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit'
            download_file(genome_url, genome_dir)
        if not os.path.isfile(chrom_len_file):
            chom_len_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
            download_file(chom_len_url, genome_dir)
    elif genome == 'mouse':
        genome_file = f'{genome_dir}mm10.2bit'
        chrom_len_file = f'{genome_dir}mm10.chrom.sizes'
        if not os.path.isfile(genome_file):
            genome_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit'
            download_file(genome_url, genome_dir)
        if not os.path.isfile(chrom_len_file):
            chom_len_url = 'https://hgdownload.soe.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes'
            download_file(chom_len_url, genome_dir)
    else:
        if not os.path.isfile(genome):
            print(f'Error: please provide 2bit genome file path')
        else:
            genome_file = genome
        if not os.path.isfile(chrom_len_file):
            print(f'Error: please provide a chromosome length file path')
        else:
            chrom_len_file = chrom_len_file
    return (genome_file, chrom_len_file)
