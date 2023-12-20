#!/usr/bin/env python3
import pandas as pd
import argparse
import ubergauss.tools as ut
from  sklearn.ensemble import RandomForestClassifier as RF
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-ih", "--human_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/human/training_data_human_context_150.npz")
    parser.add_argument("-in", "--mouse_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/mouse/training_data_mouse_context_150.npz")

    args = parser.parse_args()
    human_file = args.human_input_file
    mouse_file = args.mouse_input_file
    out_path = './../plots/features/'


    human_file = ut.nloadfile(human_file)
    X,y, ftname, _ = human_file


   # Convert the numpy array to a pandas DataFrame
    df = pd.DataFrame(X, columns=ftname)

# Calculate the correlation matrix
    corr_matrix = df.corr()

# Cluster the features
    corr_linkage = hierarchy.ward(corr_matrix)
    dendro = hierarchy.dendrogram(corr_linkage, labels=corr_matrix.columns, no_plot=True)
    # sorted_corr_index = np.array(dendro['ivl'], dtype="int")

# Create a heatmap with the clustered correlation matrix
    plt.figure(figsize=(10, 8))
    #sns.heatmap(corr_matrix.iloc[sorted_corr_index, sorted_corr_index], annot=True, fmt=".2f", cmap='coolwarm',
    #        xticklabels=corr_matrix.columns[sorted_corr_index], yticklabels=corr_matrix.columns[sorted_corr_index])
    sns.heatmap(corr_matrix)
# Show the plot
    plt.title('Clustered Feature Correlation Matrix')
    plt.savefig('feature_correlation_matrix.png', dpi=300) # Saving as PNG with high resolution







if __name__ == '__main__':
    main()
