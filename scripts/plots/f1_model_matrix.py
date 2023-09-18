import sys
import argparse
import pandas as pd
from sklearn.metrics import plot_confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import plot_roc_curve
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
import rrieval.lib as rl

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
    #print(f'AUC_PCR: {auc_prc}')
    return f1, auc_prc

def calculate_diag(name, input_path, context=150):
    df_val0 = pd.read_csv((f'{input_path}/{name}/model/{name}_context_{context}_fold0.csv'))
    df_val1 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold1.csv'))
    df_val2 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold2.csv'))
    df_val3 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold3.csv'))
    df_val4 = pd.read_csv(( f'{input_path}/{name}/model/{name}_context_{context}_fold4.csv'))

    df_cv = pd.concat([df_val0, df_val1, df_val2, df_val3, df_val4], ignore_index=True)
    # print(df_cv)

    f1, auc_prc = calculate_measures(df_cv, read=False)

    return f1, auc_prc


def generate_df(in_dict, feature_file_names):
    df_temp = pd.DataFrame(in_dict)
    # sort the colums alphabeticly
    df_sorted = df_temp.reindex(sorted(df_temp.columns), axis=1)
    # add the sorted list of dataset/model names to the datafram
    df_sorted['names']= sorted(feature_file_names)
    # make the names as index
    df_out = df_sorted.set_index('names')
    df_out.round(2)
    return df_out, df_sorted

def construct_cross_validation_call(name1,name2,file,context=150,st='on'):

    test_data = f'{file}/{name1}/feature_files/feature_filtered_{name1}_context_{context}_pos_occ'
    out = f'{file}'
    model_file = f'{file}/{name2}/model/full_{name2}_context_{context}.model'
    feat_file = f'{file}/{name2}/feature_files/training_data_{name2}_context_{context}.npz'


    call = (f'cherri eval -i1 {test_data} -g not_needed -l not_needed '
            f'-o {file} -n eval_{name2}_using_{name1} -c {context} -st {st} '
            f'-m {model_file} -mp {feat_file} -j 7 -on evaluation -ef on')

    return call


def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--input_path", required=True, help= "Path to data files")
    parser.add_argument("-st", "--structure", required=True, help= "set to on or off")
    args = parser.parse_args()

    # model type:
    input_path = args.input_path

    feature_file_names = ['human','mouse', 'human_rbp' , 'Full']
    context = 150
    st = args.structure

    #feature_file_names = ['PARIS_human','PARIS_human_RBP','PARIS_mouse','SPLASH_human','Full_human_RRIs']




#### compute the cross model evaluation!
    for name_model in feature_file_names:
        #print(f'evaluation calls for {name_model}')
        for name_test in feature_file_names:
            if name_model != name_test:
                call = construct_cross_validation_call(name_test,name_model,
                                                       input_path, context, st)
                print(f'call for {name_model} model on {name_test} data:')
                #print(call)
                rl.call_script(call)

            if name_model == name_test:
                print('implement cross validation calls later')




    # model -> [f1 socres of all datasets]
    f1_dict = {}
    AUC_dict = {}
    eval_path = f'{input_path}/evaluation/evaluation/'

    for name in feature_file_names:
        print(f'investation {name}')
        temp_dict_f1 = {}
        temp_dict_auc = {}
        for name2 in feature_file_names:
            if name != name2:
                # print(f'\n******combined with {name2}')
                key = name2 + '$' + name
                # file example: evaluation_results_PARIS_mouse_PARIS_human_RBP.cvs

                file_name = f'{eval_path}evaluation_results_eval_{name}_using_{name2}.csv'
                # print(f'eval file\n{file_name}')
                f1_cross_model, auc_prc_model = calculate_measures(file_name)
                # print(f'**done**\n')

                print(f'F1 for {name}_using_{name2}: {f1_cross_model}')

                #print(key, val[0])
                temp_dict_f1[key]= f1_cross_model
                temp_dict_auc[key]= auc_prc_model
        key_diag = name + '$' + name
        #print(f'\n******combined with {name}')
        f1, auc = calculate_diag(name, input_path, context)
        print(f'F1 for {name}: {f1}')

        # print(f'**done**\n')
        temp_dict_f1[key_diag]= f1
        temp_dict_auc[key_diag]= auc
        print(sorted(temp_dict_f1.keys()))
        f1_dict[name]=[temp_dict_f1[i] for i in sorted(temp_dict_f1.keys())]
        AUC_dict[name]=[temp_dict_auc[i] for i in sorted(temp_dict_auc.keys())]


    print(f1_dict)
    df_f1, df_sorted_f1 = generate_df(f1_dict, feature_file_names)
    print(df_sorted_f1)
    # print(df_f1)
    df_auc, df_sorted_auc = generate_df(AUC_dict, feature_file_names)
    ##df_f1.to_latex(f'{input_path}/f1_talbel', float_format="{:0.2f}".format)
    print(f'{input_path}/f1_talbel')
    df_f1.style.to_latex(f'{input_path}/f1_table')
    df_f1.to_csv(f'{input_path}/f1_table.csv',index=False)

    df_auc.style.to_latex(f'{input_path}/auc_table')
    df_auc.to_csv(f'{input_path}/auc_table.csv',index=False)

    ##for key in measures_dict:
    #    #print(key)
    #    #print(measures_dict[key][0])



if __name__ == '__main__':
    main()
