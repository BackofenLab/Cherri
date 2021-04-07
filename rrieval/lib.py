#general imports

import pandas as pd

# training imports

import csv
import pandas as pd
import pandas_profiling
import sklearn as sk

###

def read_chira_data(in_file, header='no'):
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
    #Inject labels
    pos_df['label'] = 1
    neg_df['label'] = 0
    #Dataset initial characterisation
    pos_report=pandas_profiling.ProfileReport(pos_df,title="Positive data Report")
    neg_report=pandas_profiling.ProfileReport(neg_df,title="Negative data Report")
    pos_report.to_file("positive_report.html")
    neg_report.to_file("negative_report.html")
    print(pos_df.dtypes)
    print(neg_df.dtypes)
    print(pd.get_dummies(pos_df))
    print(pd.get_dummies(neg_df))
    #Concat datasets
    ia_df = pd.concat(pos_df,neg_df)
    y = ia_df.label
    X = ia_df.drop(columns="label")
    #Create training and test dataset
    X_training, X_test, y_training, y_test = model.selection.train_test_split(X, y, test_size=0.3, random_state=42)
    #comparison dummy model
    cm = DummyClassifier()
    cm.fit(X_train, y_train)
    comparison_score = cm.score(X_test, y_test)
    print(comparison_score)


    return ""
