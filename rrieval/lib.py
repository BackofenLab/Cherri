import pandas as pd

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
