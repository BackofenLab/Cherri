#!/usr/bin/env python3
import pandas as pd
import argparse
import ubergauss.tools as ut
from  sklearn.ensemble import RandomForestClassifier as RF
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def get_feature_importances(file_zip):
    file = ut.nloadfile(file_zip)
    X,y, ftname, _ = file
    m = RF().fit(X,y)
    feature_importance = m.feature_importances_
    return feature_importance, ftname

    feature_importance, ftname = get_feature_importances()



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-ih", "--human_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/human/training_data_human_context_150.npz")
    parser.add_argument("-ir", "--human_rbp_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/human_rbp/training_data_human_rbp_context_150.npz")
    parser.add_argument("-in", "--mouse_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/mouse/training_data_mouse_context_150.npz")
    parser.add_argument("-if", "--full_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_without_graph_features/Full/training_data_Full_context_150.npz")

    args = parser.parse_args()
    # get input files
    human_file = args.human_input_file
    human_rbp_file = args.human_rbp_input_file
    mouse_file = args.mouse_input_file
    full_file = args.full_input_file

    out_path = './../plots/features/'


    # read data
    human_importance, ftname = get_feature_importances(human_file)
    human_rbp_importance, ftname = get_feature_importances(human_rbp_file)
    mouse_importance, ftname = get_feature_importances(mouse_file)
    full_importance, ftname = get_feature_importances(full_file)

    # Setting the positions and width for the bars


    # Generating names for the data points (assuming 26 points)
    names = [f'Point {chr(65+i)}' for i in range(len(human_importance))]

    # Creating a DataFrame
    df_importance_temp = pd.DataFrame({
        'Features': ftname,
        'human': human_importance,
        'mouse': mouse_importance
    })

    # Melting the DataFrame for easier plotting with seaborn
    df_importance = df_importance_temp.melt(id_vars='Features', var_name='Dataset', value_name='Value')
    print(df_importance_temp)
    print(df_importance)

    # Pivot the DataFrame to have features as indexes, datasets as columns, and values as cell values
    df_pivot = df_importance.pivot(index='Features', columns='Dataset', values='Value')

    # Calculate the difference in values between dataA and dataB
    df_pivot['difference'] = (df_pivot['human'] - df_pivot['mouse']).abs()

    # Sort the DataFrame based on the difference
    df_pivot_sorted = df_pivot.sort_values(by='difference', ascending=False)

    # Reset index to make Feature_name a column again and prepare for plotting
    df_sorted_for_plot = df_pivot_sorted.reset_index().melt(id_vars='Features', value_vars=['human', 'mouse'])

# Plot
    # Plotting
    plt.figure(figsize=(15, 8))
    sns.set_context("notebook")
    sns.set_theme(style="whitegrid")
    sns.barplot(data=df_sorted_for_plot, x='Features', y='Value', hue='Dataset')
    plt.xticks(rotation=90)
    plt.title('Compare feature importance')
    plt.xlabel('Data Points')
    plt.ylabel('Values')
    plt.tight_layout()
    plt.show()
    plt.savefig(out_path + 'feature_importance_barplot_mouse_human.png', dpi=300) # Saving as PNG with high resolution



    # Plotting
    #plt.figure(figsize=(15, 8))
    #sns.set_context("notebook")
    #sns.set_theme(style="whitegrid")
    #sns.barplot(x='Features', y='Value', hue='Dataset', data=df_importance)
    #plt.xticks(rotation=90)
    #plt.title('Compare feature importance')
    #plt.xlabel('Data Points')
    #plt.ylabel('Values')
    #plt.tight_layout()
    #plt.show()
    #plt.savefig(out_path + 'feature_importance_barplot_mouse_human.png', dpi=300) # Saving as PNG with high resolution

    #print(len(human_importance))
    #print(len(mouse_importance))
    #print(len(human_rbp_importance))
    print(len(ftname))
    print(ftname)

    df_importance_temp_all = pd.DataFrame({
        'Features': ftname,
        'human': human_importance,
        'mouse': mouse_importance,
        'human_rbp': human_rbp_importance,
        'Full': full_importance
    })

    # Melting the DataFrame for easier plotting with seaborn
    df_importance_full = df_importance_temp_all.melt(id_vars='Features', var_name='Dataset', value_name='Value')

    # Plotting
    plt.figure(figsize=(15, 8))
    sns.set_context("notebook")
    sns.barplot(x='Features', y='Value', hue='Dataset', data=df_importance_full)
    plt.xticks(rotation=90)
    plt.title('Compare feature importance')
    plt.xlabel('Data Points')
    plt.ylabel('Values')
    plt.tight_layout()
    plt.show()
    plt.savefig(out_path + 'feature_importance_barplot_full.png', dpi=300) # Saving as PNG with high resolution





    #pos = np.arange(len(human_importance))
    #bar_width = 0.35

    # Plotting the bars
    #plt.bar(pos, human_importance, bar_width, label='human')
    #plt.bar(pos + bar_width, mouse_importance, bar_width, label='mouse')

    # Adding labels and title
    #plt.xlabel('Categories')
    #plt.ylabel('Values')
    #plt.title('Side by Side Bar Plot')
    #plt.xticks(pos + bar_width / 2, ftname)
    #plt.legend()


    # Displaying the plot
    #plt.savefig(out_path + 'feature_importance_barplot_full.png', dpi=300) # Saving as PNG with high resolution






if __name__ == '__main__':
    main()
