#general imports

import pandas as pd
import subprocess
import random
import re
import os
# training imports
import csv
import pandas as pd
import pandas_profiling
import sklearn as sk
from sklearn import model_selection
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import (RandomForestClassifier)

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




def calculate_overlap(s1,e1,s2,e2):
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
    return compinde_overlap


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

        Raises
        ------
        nothing

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

        Raises
        ------
        nothing

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



#Functions for model training
def train_model(in_positive_data_filepath,in_negative_data_filepath,output_path):
    pos_df = pd.read_csv(in_positive_data_filepath, sep=',')
    neg_df = pd.read_csv(in_negative_data_filepath, sep=',')
    #Inject labels
    pos_df['label'] = 1
    neg_df['label'] = 0
    #Dataset initial characterisation
    pos_report=pandas_profiling.ProfileReport(pos_df,title="Positive data Report")
    neg_report=pandas_profiling.ProfileReport(neg_df,title="Negative data Report")
    pos_report.to_file(output_path"/positive_report.html")
    neg_report.to_file(output_path"/negative_report.html")
    #print(pos_df.dtypes)
    #print(neg_df.dtypes)
    #print(pd.get_dummies(pos_df))
    #print(pd.get_dummies(neg_df))
    #Concat datasets
    ia_df = pd.concat([pos_df,neg_df])
    y = ia_df.label
    X = ia_df.drop(columns="label")
    #Create training and test dataset
    X_training, X_test, y_training, y_test = model_selection.train_test_split(X, y, test_size=0.3, random_state=42)
    #comparison dummy model
    cm = DummyClassifier()
    cm.fit(X_training, y_training)
    dummy_comparison_score = cm.score(X_test, y_test)
    print("Dummy score:")
    print(dummy_comparison_score)
    random_forest = RandomForestClassifier(n_estimators=100, random_state=42)
    random_forest.fit(X_training, y_training)
    random_forest_comparison_score = random_forest.score(X_test, y_test)
    print("RF score:")
    print(random_forest_comparison_score)
    return ""
