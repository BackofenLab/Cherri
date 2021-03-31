import csv
import pandas as pd
import pandas_profiling

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

#Functions for model training
def train_model(in_positive_data_filepath,in_negative_data_filepath,output_path):
    pos_df = pd.read_csv(in_positive_data_filepath, sep=',')
    neg_df = pd.read_csv(in_negative_data_filepath, sep=',')
    #Dataset initial characterisation
    pos_report=pandas_profiling.ProfileReport(pos_df,title="Positive data Report")
    neg_report=pandas_profiling.ProfileReport(neg_df,title="Negative data Report")
    pos_report.to_file("positive_report.html")
    neg_report.to_file("negative_report.html")
    print(pos_df.dtypes)
    print(neg_df.dtypes)
    print(pd.get_dummies(pos_df))
    print(pd.get_dummies(neg_df))
    #print(pos_df)
    #print(neg_df)
    return ""
