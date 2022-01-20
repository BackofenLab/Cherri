import sys
import argparse
import rrieval.lib as rl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_corr_df(pos_file, neg_file):
    X, y = rl.read_pos_neg_data(pos_file, neg_file)
    df_corr = X.corr(method='spearman')
    return df_corr

def plot_clustermap(df_corr, name, out_dir):

    sns.set_theme(color_codes=True)
    ax = sns.clustermap(df_corr)

    file_name = out_dir + name + '_clustermap.pdf'
    plt.savefig(file_name)


def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-ip", "--in_positive_data_filepath", required=True, help= "Path to positive dataset")
    parser.add_argument("-in", "--in_negative_data_filepath", required=True, help= "Path to negative dataset")
    parser.add_argument("-n", "--data_name", required=True, help= "Path to negative dataset")
    parser.add_argument("-d", "--output_path", required=True, help= "Path where to store the plot")
    args = parser.parse_args()

    # model type:
    out_dir = args.output_path + '/'
    compare = True



    # data:
    if compare == True:
        #pos_human = '/vol/scratch/data/features_files/full_maxloop3/features_HEK_con150/feature_filtered_paris_HEK293T_context_150_pos_occ_pos.csv'
        #neg_human = '/vol/scratch/data/features_files/full_maxloop3/features_HEK_con150/feature_filtered_paris_HEK293T_context_150_pos_occ_neg.csv'
        #pos_mouse = '/vol/scratch/data/features_files/full_mouse_c150_maxloop3/feature_mouse_con150/feature_filtered_paris_mouse_context_150_pos_occ_pos.csv'
        #neg_mouse = '/vol/scratch/data/features_files/full_mouse_c150_maxloop3/feature_mouse_con150/feature_filtered_paris_mouse_context_150_pos_occ_neg.csv'
        #pos_humanRBP = '/vol/scratch/data/features_files/full_c150_shortRBPs/human_sortRBP/feature_filtered_paris_HEK_context_150_pos_occ_pos.csv'
        #neg_humanRBP = '/vol/scratch/data/features_files/full_c150_shortRBPs/human_sortRBP/feature_filtered_paris_HEK_context_150_pos_occ_neg.csv'

        pos_human = '/vol/scratch/data/features_files/full_maxloop3/feature_filtered_paris_HEK293T_context_150_pos_occ_pos.csv'
        neg_human = '/vol/scratch/data/features_files/full_maxloop3/feature_filtered_paris_HEK293T_context_150_pos_occ_neg.csv'
        pos_mouse = '/vol/scratch/data/features_files/full_mouse_c150_maxloop3/feature_filtered_paris_mouse_context_150_pos_occ_pos.csv'
        neg_mouse = '/vol/scratch/data/features_files/full_mouse_c150_maxloop3/feature_filtered_paris_mouse_context_150_pos_occ_neg.csv'
        pos_humanRBP = '/vol/scratch/data/features_files/full_c150_shortRBPs/feature_filtered_paris_HEK_context_150_pos_occ_pos.csv'
        neg_humanRBP = '/vol/scratch/data/features_files/full_c150_shortRBPs/feature_filtered_paris_HEK_context_150_pos_occ_neg.csv'

        df_human_corr = get_corr_df(pos_human, neg_human)
        df_mouse_corr = get_corr_df(pos_mouse, neg_mouse)
        df_humanRBP_corr = get_corr_df(pos_humanRBP, neg_humanRBP)

        plot_clustermap(df_human_corr, 'human', out_dir)
        plot_clustermap(df_mouse_corr, 'mouse', out_dir)
        plot_clustermap(df_humanRBP_corr, 'humanRBP', out_dir)

        df_diff_human = df_human_corr - df_humanRBP_corr
        df_diff_human_mouse = df_human_corr - df_mouse_corr
        #print(df_human_corr)
        #print(df_humanRBP_corr)
        #print(df_diff)

        plot_clustermap(df_diff_human, 'human_with_without_RBPs_diff', out_dir)
        plot_clustermap(df_diff_human_mouse, 'human_mouse_diff', out_dir)

        # calculate difference



    else:
        df_corr = get_corr_df(args.in_positive_data_filepath,args.in_negative_data_filepath)
        plot_clustermap(df_corr, args.data_name, out_dir)





if __name__ == '__main__':
    main()
