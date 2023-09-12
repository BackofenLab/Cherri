#!/usr/bin/env python3
import sys
import argparse
import rrieval.lib as rl
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

def delete_spaces_in_df(df):
    df_space_free = df.copy()
    df_space_free.columns = df_space_free.columns.str.strip()
    return df_space_free

#def calculate_diag(name, input_path):
#    df_val0_temp = pd.read_csv((input_path + name +'_'  + name + '_fold0.csv'))
#    df_val0 = delete_spaces_in_df(df_val0_temp)
#    df_val1_temp = pd.read_csv((input_path + name +'_' + name + '_fold1.csv'))
#    df_val1 = delete_spaces_in_df(df_val1_temp)
#    df_val2_temp = pd.read_csv((input_path + name +'_'  + name + '_fold2.csv'))
#    df_val2 = delete_spaces_in_df(df_val2_temp)
#    df_val3_temp = pd.read_csv((input_path + name +'_'  + name + '_fold3.csv'))
#    df_val3 = delete_spaces_in_df(df_val3_temp)
#    df_val4_temp = pd.read_csv((input_path + name +'_' + name + '_fold4.csv'))
#    df_val4 = delete_spaces_in_df(df_val4_temp)
#    #print(df_val0.columns)
#    #print(df_val0['true_label'])
#    #print(df_val1['true_label'])
#    #print(df_val2['true_label'])
#    #print(df_val3['true_label'])
#    #print(df_val4['true_label'])
#    df_cv = pd.concat([df_val0, df_val1, df_val2, df_val3, df_val4], ignore_index=True)
#    return df_cv



def calculate_diag(name, input_path, context=150):
    df_val0 = pd.read_csv((f'{input_path}/{name}/model/{name}_context_{context}_fold0.csv'))
    df_val1 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold1.csv'))
    df_val2 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold2.csv'))
    df_val3 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold3.csv'))
    df_val4 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold4.csv'))
    # print(df_val0.columns)
    #print(df_val0['true_label'])
    #print(df_val1['true_label'])
    #print(df_val2['true_label'])
    #print(df_val3['true_label'])
    #print(df_val4['true_label'])
    df_cv = pd.concat([df_val0, df_val1, df_val2, df_val3, df_val4], ignore_index=True)
    # print(df_cv)

    return df_cv



def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--input_path", required=True, help= "Path to data files")
    args = parser.parse_args()

    # model type:
    print('****\nRunning Precision-Recall curve calculation\n****\n')
    input_path = args.input_path

    cherri_model_path = f'{input_path}/evaluation/evaluation/'


    feature_file_names = ['human','mouse']
    feature_file_names_mfe = ['human_MFE','mouse_MFE']
    color_list = ['#F0E442','#D55E00', '#BBBBBB','#009E73','#0072B2','#808080']
    index = 0

    for name in feature_file_names:

        key_diag = name + '_' + name
        df_cv = calculate_diag(name, input_path)
        #print(df_cv.info())
        #print(df_cv['true_label'])
        precision, recall, thresholds, auc_prc = compute_prc(df_cv['true_label'].tolist(),df_cv['instance_score'].tolist())

        # for name finde E value!
        feature_path = f'{input_path}/{name}_MFE/feature_files'
        pos_feature_file = f'{feature_path}/feature_filtered_{name}_MFE_context_150_pos_occ_pos.csv'
        neg_feature_file = f'{feature_path}/feature_filtered_{name}_MFE_context_150_pos_occ_neg.csv'
        X,y = rl.read_pos_neg_data(pos_feature_file, neg_feature_file)
        precision_E, recall_E, thresholds_E, auc_prc_E = compute_prc(y.tolist(), X['E'].tolist())


        base = len(df_cv[df_cv.true_label == 1])/(len(df_cv['true_label']))

        # generat the legend information
        #label_main = 'Cherri model AUC/F1: %.2f/%.2f' % (auc_prc_main, f1_score_main)
        #label_main_seq = 'Cherri SeqStruc model AUC/F1: %.2f/%.2f' % (auc_prc_main_SeqStruc, f1_score_main_SeqStruc)
        data_name = name.split('_')[-1]
        label_Cherri = f'CheRRI {data_name} model AUC: {auc_prc.round(2)}'
        label_E = f'MFE {data_name} AUC: {auc_prc_E.round(2)}'

        #plt.plot(recall_main, precision_main, color='#009E73', label=label_main) # green
        #plt.plot(recall_main_SeqStruc, precision_main_SeqStruc, color='#0072B2', label=label_main_seq) #blue
        plt.plot(recall, precision, color=color_list[index], label=label_Cherri) #yellow
        plt.plot(recall_E, precision_E, color=color_list[(index+1)], label=label_E) #red
        #plt.plot(recall_E, precision_E, color='#E69F00', label=label_E) #orange

        plt.plot([0, 1], [base, base], color=color_list[(index+2)], linestyle='--') #grey #808080

        plt.xlabel('Recall', fontsize=14)
        plt.ylabel('Precision', fontsize=14)
        plt.xticks(fontsize=12, rotation=30)
        plt.yticks(fontsize=12)
        #plt.title('Precision Recall Characteristic (prc) Curve ' + species)
        # plt.title(input_path.split(name)[-1])

        #fontsize=20
        plt.legend(prop={'size': 12})
        #plt.show()

        plt.savefig(cherri_model_path + 'prc_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')
        index +=3
    print(f'output is stored here: {cherri_model_path}')


if __name__ == '__main__':
    main()
