#!/usr/bin/env python3
import pandas as pd
import argparse
import ubergauss.tools as ut
from  sklearn.ensemble import RandomForestClassifier as RF
import matplotlib.pyplot as plt
import numpy as np



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-ih", "--human_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_with_graph_features/human/training_data_human_context_150.npz")
    parser.add_argument("-in", "--mouse_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_with_graph_features/mouse/training_data_mouse_context_150.npz")

    args = parser.parse_args()
    human_file = args.human_input_file
    mouse_file = args.mouse_input_file
    out_path = './../plots/features/'


    human_file = ut.nloadfile(human_file)
    X,y, ftname, _ = human_file
    m = RF().fit(X,y)
    human_importance = m.feature_importances_

    imp = list(zip(human_importance, ftname))
    imp.sort()
    print(imp[-10:])
    breakpoint()
    mouse_file = ut.nloadfile(mouse_file)
    X,y, ftname, _ = mouse_file
    m = RF().fit(X,y)
    mouse_importance = m.feature_importances_

    print(human_importance)
    print(mouse_importance)

    # Setting the positions and width for the bars
    pos = np.arange(len(human_importance))
    bar_width = 0.35

    # Plotting the bars
    plt.bar(pos, human_importance, bar_width, label='human')
    plt.bar(pos + bar_width, mouse_importance, bar_width, label='mouse')

    # Adding labels and title
    plt.xlabel('Categories')
    plt.ylabel('Values')
    plt.title('Side by Side Bar Plot')
    plt.xticks(pos + bar_width / 2, ftname)
    plt.legend()

    iu

    # Displaying the plot
    plt.savefig('side_by_side_barplot.png', dpi=300) # Saving as PNG with high resolution






if __name__ == '__main__':
    main()
