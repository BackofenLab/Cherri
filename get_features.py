#!/usr/bin/env python
import pandas as pd
import math
import matplotlib as mpl
import argparse
import rrieval.lib as rl

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
    #interacting_seq = row['subseqDP']
    no_A = interacting_seq.count('A')
    no_U = interacting_seq.count('U')
    no_G = interacting_seq.count('G')
    no_C = interacting_seq.count('C')

    assert (no_A + no_U + no_G + no_C + 1) == len(interacting_seq), "something went wrong detecting the nucleotieds"
    GC_content = (no_G + no_C)/len(interacting_seq)
    return GC_content


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

    df_rri = rl.read_chira_data(input, header='yes', separater=",")
    #print(df_rri.info())
    #print(df_rri['subseqDP'])


    # MFE: mfe = df_rri['E']

    # Maximal length of the two interacting subsequence
    df_rri['len_interaction_target'] = df_rri['end1'] - df_rri['start1']
    df_rri['len_interaction_query'] = df_rri['end2'] - df_rri['start2']

    # Number of base pairs within the top 1 RRI ?
    df_rri['no_bps'] = df_rri['hybridDP'].apply(lambda x: sequence_length(x))

    # Maximal length of an interacting subsequence normalized by the number of base pairs within the top 1 RRI
    df_rri['max_inter_len'] = df_rri[['len_interaction_target', 'len_interaction_query']].max(axis=1)
    #print(df_rri['len_interaction_target'])
    #print(df_rri['len_interaction_query'])
    #print(df_rri['max_inter_len'])
    df_rri['inter_len_normby_bp'] = df_rri['max_inter_len']/df_rri['no_bps']

    # Number of base pairs within the interaction vs. the normalized maximal length of the top 1 RRI
    df_rri['bp_normby_inter_len'] = df_rri['no_bps']/df_rri['max_inter_len']

    # GC-content within interaction side
    df_rri['GC_content'] = df_rri['subseqDP'].apply(lambda x: get_GC_content(x))

    # Minimum free energy normalized by the GC-content
    df_rri['mfe_normby_GC'] = df_rri['E']/df_rri['GC_content']



    print(feature_set_list)
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
