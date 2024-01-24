#!/usr/bin/env python3
import pandas as pd
import argparse
import ubergauss.tools as ut
from  sklearn.ensemble import RandomForestClassifier as RF
import matplotlib.pyplot as plt
import numpy as np



def compute_faction(sort_list, top, report=False):
    count = 0
    #print(sort_list)
    # breakpoint()
    for i in sort_list[-top:]:
       #print(i)
       if str.isnumeric(i[1]):
           # print(f'found int')
           count += 1
    if report:
        print (f'Found {count} graph features in top {top} features (faction:{count/top})')
    # return fraction
    return count/top

def get_feature_importances(file_zip):
    file = ut.nloadfile(file_zip)
    X,y, ftname, _ = file
    m = RF(random_state=20).fit(X,y)
    feature_importance = m.feature_importances_
    return feature_importance, ftname

def compute_fraction_list(sort_list):
    fraction_list = []
    occurence = 'no'
    for pos in list(range(1, len(sort_list)+1)):
        # print(pos)
        fraction = compute_faction(sort_list, pos)
        if fraction > 0 and occurence == 'no':
            print(f'first occurence of a graph feature at position {pos}')
            occurence = 'yes'
        fraction_list.append(fraction)
    return fraction_list


def compute_fraction_for_organism(file_zip, top):
    # get feature importance
    feature_importance, ftname = get_feature_importances(file_zip)

    # sort list
    imp = list(zip(feature_importance, ftname))
    imp.sort()
    #print(imp[-top:])

    fraction = compute_faction(imp, top, True)

    fraction_list = compute_fraction_list(imp)

    return fraction_list







def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-ih", "--human_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_with_graph_features/human/training_data_human_context_150.npz")
    parser.add_argument("-in", "--mouse_input_file",
                        help= "",
                        default="/vol/scratch/data_storage/Cherri_Zenodo/Cherri_models_v4/Model_with_graph_features/mouse/training_data_mouse_context_150.npz")
    parser.add_argument("-t", "--top_feat_list",
                        help= "select number of top features",
                        default=30)

    args = parser.parse_args()
    human_file = args.human_input_file
    mouse_file = args.mouse_input_file
    top = args.top_feat_list
    out_path = './../plots/features/'

    human_fractions = compute_fraction_for_organism(human_file, top)
    print(human_fractions)

    mouse_fractions = compute_fraction_for_organism(mouse_file, top)
    print(mouse_fractions)

    positions = list(range(1,top+1))

    # Create line plots
    plt.plot(positions, human_fractions[:top], label='human')
    plt.plot(positions, mouse_fractions[:top], label='mouse')

    # Adding title and labels
    #plt.title('Line Plot with Two Lines')
    plt.xlabel('number of featrues in set')
    plt.ylabel('fraction')

    # Adding a legend to distinguish the two lines
    plt.legend()

    # Show the plot
    plt.savefig(f'{out_path}side_by_side_barplot.png', dpi=300) # Saving as PNG with high resolution




    #human_file = ut.nloadfile(human_file)
    #X,y, ftname, _ = human_file
    #m = RF().fit(X,y)
    #human_importance = m.feature_importances_

    #imp = list(zip(human_importance, ftname))
    #imp.sort()
    #print(imp[-10:])
    # breakpoit
    #mouse_file = ut.nloadfile(mouse_file)
    #X,y, ftname, _ = mouse_file
    #m = RF().fit(X,y)
    #mouse_importance = m.feature_importances_

    #print(ftname)
    #print(mouse_importance)

    #count = 0
    #for i in imp[-top:]:
    #   print(i)
    #   if str.isnumeric(i[1]):
    #       print(f'found int')
    #       count += 1
    #print (f'Found {count} graph features in top {top} features (faction:{count/top})')

    # Setting the positions and width for the bars
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
    #plt.savefig('side_by_side_barplot.png', dpi=300) # Saving as PNG with high resolution






if __name__ == '__main__':
    main()
