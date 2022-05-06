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

def get_model_pathes(model_type, model_path):
    # model_selection
    model_dummy = model_path + '/HEK_'+ model_type + '_dc.obj'
    model_svc = model_path + '/HEK_' + model_type + '_svc.obj'
    model_clf = model_path + '/HEK_' + model_type + '_clf.obj'
    return model_dummy, model_svc, model_clf


def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-ip", "--in_positive_data_filepath", required=True, help= "Path to positive dataset")
    parser.add_argument("-in", "--in_negative_data_filepath", required=True, help= "Path to negative dataset")
    parser.add_argument("-mt", "--model_type", required=True, help= "model type: E, allEs or full")
    parser.add_argument("-d", "--output_path", required=True, help= "Path where to store the plot")
    args = parser.parse_args()

    # model type:
    model_type = args.model_type
    out_dir = args.output_path + '/' + model_type + '_'
    model_path = '/vol/scratch/data_storage/output/Base_model/'
    model_name_E = 'E'
    model_name_allEs = 'allEs'
    model_name_full = 'full'


    model_dummy, model_svc, model_clf = get_model_pathes(model_type, model_path)
    X, y = rl.read_pos_neg_data(args.in_positive_data_filepath,args.in_negative_data_filepath)

    if model_type == model_name_E:
        X = X[['E']]
    elif model_type == model_name_allEs:
        X = X[['E','E_hybrid', 'ED1', 'ED2']]
    elif model_type == model_name_full:
        X = X.drop(['subseqDP', 'hybridDP'], axis=1)
    else:
        print('incorrect model type only use: E, allEs or full')


    model_dummy, y_pred_dummy, score_dummy = rl.classify2(X,model_dummy, 'proba')
    model_svc, y_pred_svc, score_svc = rl.classify2(X,model_svc, 'decision')
    model_clf, y_pred_clf, score_clf = rl.classify2(X,model_clf, 'decision')




    ####Confusion matrix########################################
    #fig, (ax1, ax2) = plt.subplots(1, 2)

    #ax1.plot_confusion_matrix(estimator=model1, X=X, y_true=y,cmap='Blues')
    #ax2.plot_confusion_matrix(estimator=model2, X=X, y_true=y,cmap='Blues')

    #plt.savefig(args.output_path + 'cm.png')


    ###PRC!#################################

    df_main_model = rl.read_chira_data('/vol/scratch/data_storage/Model_with_graph_feat/evaluation/evaluation/evaluation_results_PARIS_human_PARIS_mouse.cvs', 'yes', ",")
    df_main_model_SeqStruc = rl.read_chira_data('/vol/scratch/data_storage/Model_with_graph_feat/evaluation/evaluation/evaluation_results_PARIS_human_PARIS_human_RBP.cvs', 'yes', ",")
    #print(df_main_model.info())
    #print(df_main_model.columns)
    #print(df_main_model.instance_score)
    #print(df_main_model['instance_score'].tolist())
    #print(y.tolist())


    f1_score_main = f1_score(df_main_model['true_label'].tolist(), df_main_model['predicted_label'].tolist())
    f1_score_main_SeqStruc = f1_score(df_main_model_SeqStruc['true_label'].tolist(), df_main_model_SeqStruc['predicted_label'].tolist())
    f1_score_dv = f1_score(y.tolist(), y_pred_dummy)
    f1_score_svc = f1_score(y.tolist(), y_pred_svc)
    f1_score_lr = f1_score(y.tolist(), y_pred_clf)

    #precision_main, recall_main, thresholds_main, auc_prc_main = compute_prc(df_main_model['true_label'].tolist(), df_main_model['instance_score'].tolist())
    #precision_main_SeqStruc, recall_main_SeqStruc, thresholds_main_SeqStruc, auc_prc_main_SeqStruc = compute_prc(df_main_model_SeqStruc['true_label'].tolist(), df_main_model_SeqStruc['instance_score'].tolist())
    precision_svc, recall_svc, thresholds_svc, auc_prc_svc = compute_prc(y.tolist(), score_svc)
    precision_lr, recall_lr, thresholds_lr, auc_prc_lr = compute_prc(y.tolist(), score_clf)
    precision_E, recall_E, thresholds_E, auc_prc_E = compute_prc(y.tolist(), X['E'].tolist() )


    base = len(df_main_model[df_main_model.true_label == 1])/(len(df_main_model['true_label']))

    # generat the legend information
    #label_main = 'Cherri model AUC/F1: %.2f/%.2f' % (auc_prc_main, f1_score_main)
    #label_main_seq = 'Cherri SeqStruc model AUC/F1: %.2f/%.2f' % (auc_prc_main_SeqStruc, f1_score_main_SeqStruc)
    label_scv = 'Base SVC model AUC/F1: %.2f/%.2f' % (auc_prc_svc, f1_score_svc)
    label_lr = 'Base LR model AUC/F1: %.2f/%.2f' % (auc_prc_lr,f1_score_lr)
    label_E = 'MFE AUC: %.2f' % (auc_prc_E)

    #plt.plot(recall_main, precision_main, color='#009E73', label=label_main) # green
    #plt.plot(recall_main_SeqStruc, precision_main_SeqStruc, color='#0072B2', label=label_main_seq) #blue
    plt.plot(recall_svc, precision_svc, color='#F0E442', label=label_scv) #yellow
    plt.plot(recall_lr, precision_lr, color='#D55E00', label=label_lr) #red
    plt.plot(recall_E, precision_E, color='#E69F00', label=label_E) #orange

    plt.plot([0, 1], [base, base], color='#BBBBBB', linestyle='--') #grey #808080

    plt.xlabel('Recall', fontsize=14)
    plt.ylabel('Precision', fontsize=14)
    plt.xticks(fontsize=12, rotation=30)
    plt.yticks(fontsize=12)
    #plt.title('Precision Recall Characteristic (prc) Curve ' + species)
    plt.title(model_type)

    #fontsize=20
    plt.legend(prop={'size': 12})
    #plt.show()

    plt.savefig(out_dir+'prc_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')


    ####################################################################


    #model1_disp = plot_roc_curve(model1, X, y)
    #model2_disp = plot_roc_curve(model2, X, y, ax=model1_disp.ax_)
    #model3_disp = plot_roc_curve(model3, X, y, ax=model2_disp.ax_)
    #model1_disp.figure_.suptitle("ROC curve comparison")
    #plt.savefig(args.output_path + 'roc.png')


    # compute feature importance!


if __name__ == '__main__':
    main()
