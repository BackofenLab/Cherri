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
import pandas as pd
import numpy as np
import pandas_profiling
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
#import xgboost
import pickle
from sklearn.linear_model import Lasso
#import seaborn as sns
###

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
        df_temp = pd.read_table(in_file, header=None, sep=separater)
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
    elif header == 'yes':
        df_interactions = pd.read_table(in_file, sep=separater)
    return df_interactions


def call_script(call,reprot_stdout=False):
    """
    starts a subprosses to call a script and checks for errors.


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
    out, err = process.communicate()
    #print(out.decode('utf-8'))
    error = err.decode('utf-8')

    assert not error, "script is complaining:\n%s\n%s" %(call, error)
    if reprot_stdout == True:
        # out = out.decode('utf-8')
        return out

def calculate_overlap(s1,e1,s2,e2,len_flag=False):
    """
    Building for each replicat a inter object

        Parameters
        ----------
        s1: start of one sequence of the first replicat
        e1: end of one sequence of the first replicat
        s2: start of one sequence of the current replicat
        e2: end of one sequence of the current replicat


        Raises
        ------
        nothing

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



def get_chrom_list_no_numbers(df_interactions, chrom):
    """
    Generates a unique list of chromosmes or conticts for both interaction
    partners

        Parameters
        ----------
        df_interactions : df including the filtered RRIs
        chrom :  string indicating from wich seq the chromosome is


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
        seqs_dic = read_fasta_into_dic(out_fa, skip_n_seqs=True)
        # Output RNA sequences.
        #fasta_output_dic(seqs_dic, out_fa, split=True)
    return seqs_dic


def check_context(df, seq_tag, chrom_dict):
    """
    check that the extende contxt is not to short or long!

        Parameters
        ----------
        df: bed df


        Returns
        -------
        df
            df with changed postions
    chrom_dict = {'chr1':60,
    ...         'chr2':80}
    >>> data = {'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'chr1', 'chr1', 'chr1'],
    ...         'end_1st':[30, 60, 70, 70],
    ...         'start_2end':[0, -40, 40, -2],
    ...         'end_2end':[30, 60, 70, 70],
    ...         'chrom_2end':['chr1', 'chr1', 'chr1', 'chr1']}
    >>> df = pd.DataFrame(data)
    >>> check_context(df, 'target', chrom_dict)
    pd.DataFrame({'start_1st':[0, 0, 40, 0],
    ...         'chrom_1st':['chr1', 'chr1', 'chr1', 'chr1'],
    ...         'end_1st':[30, 60, 60, 60],
    ...         'start_2end':[0, -40, 40, -2],
    ...         'end_2end':[30, 60, 70, 70],
    ...         'chrom_2end':['chr1', 'chr1', 'chr1', 'chr1']})
    >>> check_context(df, 'query', chrom_dict)
    pd.DataFrame({'start_1st':[0, -40, 40, -2],
    ...         'chrom_1st':['chr1', 'chr1', 'chr1', 'chr1'],
    ...         'end_1st':[30, 60, 70, 70],
    ...         'start_2end':[0, 0, 40, 0],
    ...         'end_2end':[30, 60, 60, 60],
    ...         'chrom_2end':['chr1', 'chr1', 'chr1', 'chr1']})
        """
    no_seq_out_boder = 0
    if seq_tag == 'target':
        no_seq_out_boder += len(df[df.start_1st <=0])
        no_seq_out_boder += len(df[df.end_1st >= df['chrom_1st'].apply(lambda x: chrom_dict[x])])
        df.loc[df.start_1st <= 0, 'start_1st'] = 0
        df.loc[df.end_1st >= df['chrom_1st'].apply(lambda x: chrom_dict[x]), 'end_1st'] = df['chrom_1st'].apply(lambda x: chrom_dict[x])

    elif seq_tag == 'query':
        no_seq_out_boder += len(df[df.start_2end <=0])
        no_seq_out_boder += len(df[df.end_2end >= df['chrom_2end'].apply(lambda x: chrom_dict[x])])
        df.loc[df.start_2end <= 0, 'start_2end'] = 0
        df.loc[df.end_2end >= df['chrom_2end'].apply(lambda x: chrom_dict[x]), 'end_2end'] = df['chrom_2end'].apply(lambda x: chrom_dict[x])
    print('Warning: added context to %s is out of bourder for %i instances'%(seq_tag,no_seq_out_boder))
    return df




def get_context(seq_tag, df, out_dir, in_2bit_file, context, chrom_len_file):
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
    no_del_entys = 0
    chrom_dict = read_table_into_dic(chrom_len_file)
    #print(chrom_dict)
    if seq_tag == 'target':
        df_bed = df[['chrom_1st', 'start_1st', 'end_1st', 'ID1', 'interaction_no', 'strand_1st']].copy()
        #print(df_bed.tail())
        df_bed['chrom_1st'] = df_bed['chrom_1st'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_1st', 'end_1st')
        df_context = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_target'
        col_id = 'ID1'
        df_context_filted = df_context[df_context.chrom_1st != False]
        no_del_entys += len(df_context) - len(df_context_filted)
    elif seq_tag == 'query':
        df_bed = df[['chrom_2end', 'start_2end', 'end_2end', 'ID2', 'interaction_no', 'strand_2end']].copy()
        df_bed['chrom_2end'] = df_bed['chrom_2end'].apply(lambda x: check_convert_chr_id(x))
        df_context =  add_context(df_bed, context, 'start_2end', 'end_2end')
        df_context = check_context(df_context, seq_tag, chrom_dict)
        col_name = 'con_query'
        col_id = 'ID2'
        df_context_filted = df_context[df_context.chrom_2end != False]
        no_del_entys += len(df_context) - len(df_context_filted)
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


#Functions negative data
def shuffle_sequence(seq, times, kind_of_shuffel):
    """
    shuffle on given sequence x times

        Parameters
        ----------
        seq: sequence
        times: amount of shuffling
        kind_of_shuffel: 1 -> Mononucleotide; 2 -> Dinucleotide

        Raises
        ------
        nothing

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


def bp_suffeling(hybrid_seq, IntaRNA_prediction,times):
    """
    basepair shufelling of the given IntaRNA prediction

        Parameters
        ----------
        hybrid_seq: tartet sequence
        IntaRNA_prediction: query sequence

        Raises
        ------
        nothing

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

        Raises
        ------
        nothing

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



#Functions for model training
def train_model(in_positive_data_filepath,in_negative_data_filepath,output_path):
    #pos_df = pd.read_csv(in_positive_data_filepath, sep=',')
    #neg_df = pd.read_csv(in_negative_data_filepath, sep=',')
    #Inject labels
    #pos_df['label'] = 1
    #neg_df['label'] = 0
    #Dataset initial characterisation
    #reporting=0
    #if(reporting):
        #pos_report=pandas_profiling.ProfileReport(pos_df,title="Positive data Report")
        #neg_report=pandas_profiling.ProfileReport(neg_df,title="Negative data Report")
        #pos_report.to_file(output_path + "/positive_report.html")
        #neg_report.to_file(output_path + "/negative_report.html")
    #print(pos_df.dtypes)
    #print(neg_df.dtypes)
    #print(pd.get_dummies(pos_df))
    #print(pd.get_dummies(neg_df))
    #Concat datasets
    #ia_df = pd.concat([pos_df,neg_df])
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

def classify(in_data_filepath,in_model_filepath,output_path):
    X = pd.read_csv(in_data_filepath, sep=',')
    model_handle = open(in_model_filepath,'rb')
    model = pickle.load(model_handle)
    model_handle.close()
    y_pred=model.predict(X)
    print(y_pred)
    return y_pred

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

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))

    return accuracy
