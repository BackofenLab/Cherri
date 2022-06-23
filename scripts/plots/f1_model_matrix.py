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
    f1 = f1_score(df_model['true_label'].tolist(), df_model['predicted_label'].tolist())
    #print(f1)
    precision, recall, thresholds, auc_prc = compute_prc(df_model['true_label'].tolist(), df_model['instance_score'].tolist())
    print(f'AUC_PCR: {auc_prc}')
    return f1, auc_prc

def calculate_diag(name, input_path):
    df_val0 = pd.read_csv((input_path + name +'_'  + name + '_fold0.csv'))
    df_val1 = pd.read_csv((input_path + name +'_' + name + '_fold1.csv'))
    df_val2 = pd.read_csv((input_path + name +'_'  + name + '_fold2.csv'))
    df_val3 = pd.read_csv((input_path + name +'_'  + name + '_fold3.csv'))
    df_val4 = pd.read_csv((input_path + name +'_' + name + '_fold4.csv'))
    df_cv = pd.concat([df_val0, df_val1, df_val2, df_val3, df_val4], ignore_index=True)
    f1, auc_prc = calculate_measures(df_cv, read=False)

    return f1, auc_prc


def generate_df(in_dict, feature_file_names):
    df_temp = pd.DataFrame(in_dict)
    df_sorted = df_temp.reindex(sorted(df_temp.columns), axis=1)
    df_sorted['names']= sorted(feature_file_names)
    df_out = df_sorted.set_index('names')
    df_out.round(2)
    return df_out




def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--input_path", required=True, help= "Path to data files")
    args = parser.parse_args()

    # model type:
    input_path = args.input_path


    feature_file_names = ['PARIS_human','PARIS_mouse', 'PARIS_human_RBP','Full']

    #feature_file_names = ['PARIS_human','PARIS_human_RBP','PARIS_mouse','SPLASH_human','Full_human_RRIs']



    # contains all measurments
    measures_dict_f1 = {}
    measures_dict_auc = {}
    # model -> [f1 socres of all datasets]
    f1_dict = {}
    AUC_dict = {}

    for name in feature_file_names:
        temp_dict_f1 = {}
        temp_dict_auc = {}
        for name2 in feature_file_names:
            if name != name2:
                key = name + '$' + name2
                # file example: evaluation_results_PARIS_mouse_PARIS_human_RBP.cvs
                file_name = f'{input_path}evaluation_results_{name}_{name2}.csv'
                f1_cross_model, auc_prc_model = calculate_measures(file_name)
                #measures_dict_f1[key] = f1_cross_model
                #measures_dict_auc[key] = auc_prc_model
                #print(key, val[0])
                temp_dict_f1[name2]= f1_cross_model
                temp_dict_auc[name2]= auc_prc_model
        key_diag = name + '$' + name
        f1, auc = calculate_diag(name, input_path)
        #measures_dict_f1[key_diag] = f1
        #measures_dict_auc[key_diag] = auc
        temp_dict_f1[name]= f1
        temp_dict_auc[name]= auc
        print(key_diag, f1)
        f1_dict[name]=[temp_dict_f1[i] for i in sorted(temp_dict_f1.keys())]
        AUC_dict[name]=[temp_dict_auc[i] for i in sorted(temp_dict_auc.keys())]



    df_f1 = generate_df(f1_dict, feature_file_names)
    df_auc = generate_df(AUC_dict, feature_file_names)
    #df_f1.to_latex(f'{input_path}/f1_talbel', float_format="{:0.2f}".format)
    print(f'{input_path}/f1_talbel')
    df_f1.style.to_latex(f'{input_path}/f1_table')

    df_auc.style.to_latex(f'{input_path}/auc_table')

    #for key in measures_dict:
        #print(key)
        #print(measures_dict[key][0])



if __name__ == '__main__':
    main()
