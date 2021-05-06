#!/usr/bin/env python
import pandas as pd
import math
import matplotlib as mpl
import argparse
import rrieval.lib as rl
import re
import subprocess



def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True
                                           , help= "path to input files")
    parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", required=True,
                        help= "path to output dir")
    parser.add_argument("-e", "--experiment", action="store", dest="experiment", required=True,
                         help= "prefix indicating the experiment and subest to be analyzed")



    args = parser.parse_args()

    input_dir = args.input_dir
    out_dir = args.out_dir
    experiment = args.experiment


    input_dir = input_dir + '/'
    # infos in the overview table are used to generate the calls for the model
    file_overview = input_dir + 'overview.tabular'


    df_overview = pd.read_table(file_overview)
    #print(df_overview.info())

    models = ['DummyClassifier_auc', 'LogisticRegression_auc', 'DecisionTreeClassifier_auc',
              'KNeighborsClassifier_auc', 'GaussianNB_auc','SVC_auc',
              'RandomForestClassifier_auc', 'XGBClassifier_auc']
    models_std = ['DummyClassifier_std', 'LogisticRegression_std', 'DecisionTreeClassifier_std',
              'KNeighborsClassifier_std', 'GaussianNB_std','SVC_std',
               'RandomForestClassifier_std', 'XGBClassifier_std']

    COLUMN_NAMES= list(df_overview.columns) + models + models_std
    #print(COLUMN_NAMES)
    df_result = pd.DataFrame(columns=COLUMN_NAMES)
    #print(df_result)


    for index, row in df_overview.iterrows():
        if re.match(r"[0-9]+a", row['id']):
            pos = row['id']
            continue
        modle_values = models
        models_std_val = models_std
        neg = row['id']
        call = 'python training.py -i ' + input_dir + pos + ' -n ' + input_dir + neg +' -d ' + out_dir
        print(call)
        process = subprocess.Popen(call, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, shell=True)
        for idx, line in enumerate(process.stdout):
            line = line.decode("utf-8").strip().split('\t')
            #print(line[2])
            if line[0].strip() == 'DummyClassifier':
                modle_values[0] = line[2]
                models_std_val[0] = line[4]
            elif line[0].strip() == 'LogisticRegression':
                modle_values[1] = line[2]
                models_std_val[1] = line[4]
            elif line[0].strip() == 'DecisionTreeClassifier':
                modle_values[2] = line[2]
                models_std_val[2] = line[4]
            elif line[0].strip() == 'KNeighborsClassifier':
                modle_values[3] = line[2]
                models_std_val[3] = line[4]
            elif line[0].strip() == 'GaussianNB':
                modle_values[4] = line[2]
                models_std_val[4] = line[4]
            elif line[0].strip() == 'SVC':
                modle_values[5] = line[2]
                models_std_val[5] = line[4]
            elif line[0].strip() == 'RandomForestClassifier':
                modle_values[6] = line[2]
                models_std_val[6] = line[4]
            elif line[0].strip() == 'XGBClassifier':
                modle_values[7] = line[2]
                models_std_val[7] = line[4]
            else:
                print('unknown line')

        # append values
        print(modle_values)
        #print(list(row))
        values = list(row) + modle_values + models_std_val
        print(values)
        df_result.loc[len(df_result)] = values
        print(df_result.loc[len(df_result)-1])
    print(df_result)
    df_result.to_csv(out_dir + experiment + '_modle_AUC.csv', index=False)

if __name__ == '__main__':
    main()
