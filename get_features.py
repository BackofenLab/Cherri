#!/usr/bin/env python
import pandas as pd
import math
import matplotlib as mpl
import argparse
from .rrieval import lib as rl
import re


##### possible features to add:
# -> sequence based features:
# - contents and skews https://en.wikipedia.org/wiki/GC_skew
# - https://github.com/caballero/SeqComplex
# -  complexity of the sequences: Wootton–Federhen complexity, entropy of the symbols,
#    entropy of symbols in a Markov model, linguistic complexity:
#    doi: 10.1093/nar/gku356
# -> structure based features:
# - MFE, frequency of MFE, ensemble diversity (ED), thermodynamic z-score and its p-value
#  DOI: 10.1038/s41598-017-17510-y



def sequence_length(hybrid):
    """
    compute number of base pairs

        Parameters
        ----------
        hybrid: interaction in a dot bracked notation

        Raises
        ------
        nothing

        Returns
        -------
        no_bps
            number of base pairs


        """
    hybrid = str(hybrid)
    no_bps = hybrid.count('(')
    no_bps_test = hybrid.count(')')
    #print(hybrid)
    #print(no_bps)
    #print(no_bps_test)
    #assert no_bps == no_bps_test, "not equal number of open and closing brackes: %s"%(hybrid)

    return no_bps

def count_number_of_seeds(seed_pos):
    """
    compute number of base pairs

        Parameters
        ----------
        seed_pos_target : start or end postion of tarted
        seed_pos_query : start or end postion of query

        Raises
        ------
        nothing

        Returns
        -------
        no_seed
            number of seeds


        """
    #seed_pos_target = '1:5:8:10'
    #seed_pos_query = '2:4:6:8'
    #print(seed_pos)
    match = re.search(r":", str(seed_pos))
    if match:
        #print('found pattern')
        no_seed = seed_pos.count(':') + 1
    else:
        no_seed = 1
    return no_seed


def calc_seq_entropy(seq_l, ntc_dic):
    """
    Given a dictionary of nucleotide counts for a sequence ntc_dic and
    the length of the sequence seq_l, compute the Shannon entropy of
    the sequence.

    Formula (see CE formula) taken from:
    https://www.ncbi.nlm.nih.gov/pubmed/15215465

    >>> seq_l = 8
    >>> ntc_dic = {'A': 8, 'C': 0, 'G': 0, 'U': 0}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    0
    >>> ntc_dic = {'A': 4, 'C': 4, 'G': 0, 'U': 0}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    0.5
    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'U': 2}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    1.0

    """
    # For DNA or RNA, k = 4.
    k = 4
    # Shannon entropy.
    ce = 0
    for nt in ntc_dic:
        c = ntc_dic[nt]
        if c != 0:
            ce += (c/seq_l) * math.log((c/seq_l), k)
    if ce == 0:
        return 0
    else:
        return -1*ce



def seq_count_nt_freqs(seq,
                       rna=True,
                       count_dic=False):
    """
    Count nucleotide (character) frequencies in given sequence seq.
    Return count_dic with frequencies.
    If count_dic is given, add count to count_dic.

    rna:
    Instead of DNA dictionary, use RNA dictionary (A,C,G,U) for counting.

    count_dic:
    Supply a custom dictionary for counting only characters in
    this dictionary + adding counts to this dictionary.

    >>> seq = 'AAAACCCGGT'
    >>> seq_count_nt_freqs(seq)
    {'A': 4, 'C': 3, 'G': 2, 'T': 1}
    >>> seq = 'acgtacgt'
    >>> seq_count_nt_freqs(seq)
    {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    """

    assert seq, "given sequence string seq empty"
    if not count_dic:
        count_dic = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        if rna:
            count_dic = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    # Conver to list.
    seq_list = list(seq)
    for nt in seq_list:
        if nt in count_dic:
            count_dic[nt] += 1
    return count_dic

def comput_complexity(seq):
    """
    compute number of base pairs

        Parameters
        ----------
        seq: RNA sequence

        Raises
        ------
        nothing

        Returns
        -------
        GC_content
            GC_content


        """
    seq_len = len(seq)
    count_dic = seq_count_nt_freqs(seq)
    complexity = calc_seq_entropy(seq_len, count_dic)

    return complexity


def get_GC_content(interacting_seq):
    """
    compute number of base pairs

        Parameters
        ----------
        interacting_seq: interaction sequence

        Raises
        ------
        nothing

        Returns
        -------
        GC_content
            GC_content


        """
    #print(interacting_seq)
    interacting_seq = str(interacting_seq)
    no_A = interacting_seq.count('A')
    no_U = interacting_seq.count('U')
    no_G = interacting_seq.count('G')
    no_C = interacting_seq.count('C')
    #print(no_A)

    assert (no_A + no_U + no_G + no_C + 1) == len(interacting_seq), "something went wrong detecting the nucleotieds"
    GC_content = (no_G + no_C)/len(interacting_seq)
    return GC_content

def loadDF(input, validation=0, shuffled_flag=False):
    df_rri = rl.read_chira_data(input, header='yes', separater=",")

    if validation == 1:
        df_rri['subseqDP'] = df_rri['subseq1'] + '&' + df_rri['subseq2']
        df_rri['target'] = df_rri['seq1']
        df_rri['query'] = df_rri['seq2']
    print(df_rri.info())
    #print(df_rri['subseqDP'])


    # MFE: mfe = df_rri['E']

    # Maximal length of the two interacting subsequence
    df_rri['len_interaction_target'] = df_rri['end1'] - df_rri['start1'] + 1
    df_rri['len_interaction_query'] = df_rri['end2'] - df_rri['start2'] + 1

    # Number of base pairs within the top 1 RRI ?
    if validation == 0:
        df_rri['no_bps'] = df_rri['hybridDP'].apply(lambda x: sequence_length(x))
    elif validation == 1:
        df_rri['no_bps'] = df_rri['hybridDPfull'].apply(lambda x: sequence_length(x))

    # Maximal length of an interacting subsequence normalized by the number of base pairs within the top 1 RRI
    df_rri['max_inter_len'] = df_rri[['len_interaction_target', 'len_interaction_query']].max(axis=1)
    #print(df_rri['len_interaction_target'])
    #print(df_rri['len_interaction_query'])
    #print(df_rri['max_inter_len'])
    df_rri['inter_len_normby_bp'] = df_rri['max_inter_len']/df_rri['no_bps']

    # Number of base pairs within the interaction vs. the normalized maximal length of the top 1 RRI
    df_rri['bp_normby_inter_len'] = df_rri['no_bps']/df_rri['max_inter_len']

    # GC-content within interaction side
    #print(df_rri['subseqDP'].type)

    df_rri['GC_content'] = df_rri['subseqDP'].apply(lambda x: get_GC_content(x))

    # Number of seeds seedStart1
    df_rri['no_seeds'] = df_rri['seedStart1'].apply(lambda x: count_number_of_seeds(x))

    # E_hybrid,ED1,ED2'
    df_rri['max_ED'] = df_rri[['ED1', 'ED2']].max(axis=1)

    # Energys normalized by the GC-content
    df_rri['mfe_normby_GC'] = df_rri['E']/df_rri['GC_content']
    df_rri['max_ED_normby_GC'] = df_rri['max_ED']/df_rri['GC_content']
    df_rri['E_hybrid_normby_GC'] = df_rri['E_hybrid']/df_rri['GC_content']

    # sequence complexety shannon entropy
    if shuffled_flag:
        df_rri['complex_target'] = df_rri['target'].apply(lambda x: comput_complexity(x))
        #print(df_rri['complex_target'])
        df_rri['complex_query'] = df_rri['query'].apply(lambda x: comput_complexity(x))
    else:
        # sequence complexety shannon entropy
        df_rri['complex_target'] = df_rri['con_query'].apply(lambda x: comput_complexity(x))
        #print(df_rri['complex_target'])
        df_rri['complex_query'] = df_rri['con_query'].apply(lambda x: comput_complexity(x))
        df_rri['side_target'] = df_rri['subseqDP'].apply(lambda x: x.split('&')[0])
        #print(df_rri['complex_target'])
        df_rri['side_query'] = df_rri['subseqDP'].apply(lambda x: x.split('&')[1])
        df_rri['complex_target_side'] = df_rri['side_target'].apply(lambda x: comput_complexity(x))
        #print(df_rri['complex_target'])
        df_rri['complex_query_side'] = df_rri['side_query'].apply(lambda x: comput_complexity(x))
    #print(df_rri['complex_query'])
    return df_rri





def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", action="store", dest="input", required=True
                                           , help= "path to input file")
    parser.add_argument("-f", "--feature_set_list", action="store", nargs='+',
                        dest="feature_set_list", required=True,
                        help= "set of features the script will output ")
    parser.add_argument("-o", "--output_file", action="store", dest="output_file", required=True
                                           , help= "output file path inclusive of the file name")

    args = parser.parse_args()

    input = args.input
    feature_set_list = args.feature_set_list
    output_file = args.output_file
    shuffled_flag = False
    validation = 0
    df_rri= loadDF(input, validation, shuffled_flag)
    if feature_set_list[0] == 'all' or feature_set_list[0] == 'All':
        # ToDo: now we have duplicats. find a more recine set of all!!
        df_feature = df_rri
        print('OUTput contians all featurs')
    else:
        for col_name in feature_set_list:
            print(col_name)
            #print(df_rri[col_name])
            print('#################')
        # generate output:
        df_feature = df_rri[feature_set_list].copy()
    df_feature.to_csv(output_file, index=False)

    # Go over list of featurs for the output dir and generate table:



        # Number of possible seeds



if __name__ == '__main__':
    main()
