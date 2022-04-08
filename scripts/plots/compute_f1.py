#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score


def compute_prc(labels, scores):
    precision, recall, thresholds = precision_recall_curve(
    labels, scores)
    auc_prc =metrics.auc(recall, precision)
    #auc_prc = 0
    return precision, recall, thresholds, auc_prc

def calculate_measures(df_model):
    df_model.columns = df_model.columns.str.strip()
    #print(df_model.columns)
    #print(df_model['true_label'])
    f1 = f1_score(df_model['true_label'].tolist(), df_model['prediction'].tolist())
    #print(f1)
    #precision, recall, thresholds, auc_prc = compute_prc(df_model['true_label'].tolist(), df_model['instance_score'].tolist())
    return f1

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--in_data_filepath", required=True, help= "Path to positive dataset")
    args = parser.parse_args()

    # model type:
    df_data = pd.read_csv(args.in_data_filepath)
    f1 = calculate_measures(df_data)
    print(f'F1-Score: {f1}')


if __name__ == '__main__':
    main()
