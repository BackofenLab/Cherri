import sys
import argparse
import pandas as pd
from sklearn.metrics import plot_confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import plot_roc_curve
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score


def compute_prc(labels, scores):
    precision, recall, thresholds = metrics.precision_recall_curve(
    labels, scores)
    auc_prc =metrics.auc(recall, precision)
    #auc_prc = 0
    return precision, recall, thresholds, auc_prc


def calculate_measures(data, read=True):
    if read:
        df_model = pd.read_csv(data)
    else:
        df_model = data
    df_model.columns = df_model.columns.str.strip()
    #print(df_model.columns)
    #print(df_model['true_label'])
    f1 = f1_score(df_model['true_label'].tolist(), df_model['predicted_label'].tolist(), average='macro')
    #print(f1)
    #precision, recall, thresholds, auc_prc = compute_prc(df_model['true_label'].tolist(), df_model['instance_score'].tolist())
    return f1

def calculate_diag(name, input_path):
    df_val0 = pd.read_csv((input_path + name +'_'  + name + '_fold0.csv'))
    df_val1 = pd.read_csv((input_path + name +'_' + name + '_fold1.csv'))
    df_val2 = pd.read_csv((input_path + name +'_'  + name + '_fold2.csv'))
    df_val3 = pd.read_csv((input_path + name +'_'  + name + '_fold3.csv'))
    df_val4 = pd.read_csv((input_path + name +'_' + name + '_fold4.csv'))
    df_cv = pd.concat([df_val0, df_val1, df_val2, df_val3, df_val4], ignore_index=True)
    f1 = calculate_measures(df_cv, read=False)

    return f1

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--input_path", required=True, help= "Path to data files")
    args = parser.parse_args()

    # model type:
    input_path = args.input_path


    feature_file_names = ['PARIS_human','PARIS_mouse', 'PARIS_human_RBP','Full']

    #feature_file_names = ['PARIS_human','PARIS_human_RBP','PARIS_mouse','SPLASH_human','Full_human_RRIs']



    # contains all measurments
    measures_dict = {}
    # model -> [f1 socres of all datasets]
    f1_dict = {}

    for name in feature_file_names:
        temp_dict = {}
        for name2 in feature_file_names:
            if name != name2:
                key = name + '$' + name2
                # file example: evaluation_results_PARIS_mouse_PARIS_human_RBP.cvs
                file_name = f'{input_path}evaluation_results_{name}_{name2}.cvs'
                f1_cross_model = calculate_measures(file_name)
                measures_dict[key] = f1_cross_model
                #print(key, val[0])
                temp_dict[name2]= measures_dict[key]
        key_diag = name + '$' + name
        f1 = calculate_diag(name, input_path)
        measures_dict[key_diag] = f1
        temp_dict[name]= measures_dict[key_diag]
        print(key_diag, f1)
        f1_dict[name]=[temp_dict[i] for i in sorted(temp_dict.keys())]

    # generate table form f1 dict
    df_temp = pd.DataFrame(f1_dict)
    df_sorted = df_temp.reindex(sorted(df_temp.columns), axis=1)
    df_sorted['names']= sorted(feature_file_names)
    df_f1 = df_sorted.set_index('names')
    df_f1.round(2)
    #df_f1.to_latex(f'{input_path}/f1_talbel', float_format="{:0.2f}".format)
    print(f'{input_path}/f1_talbel')
    df_f1.style.to_latex(f'{input_path}/f1_talbel')


    #for key in measures_dict:
        #print(key)
        #print(measures_dict[key][0])



if __name__ == '__main__':
    main()
