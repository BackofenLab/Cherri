#!/usr/bin/env python
import sys
import argparse
import rrieval.lib as rl
import sklearn.datasets
import sklearn.metrics
import pickle

import autosklearn.classification

#from autosklearn.experimental.askl2 import AutoSklearn2Classifier as ASK2
#from autosklearn.classification import AutoSklearnClassifier as ASK1
#import autosklearn.metrics

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--in_positive_data_filepath", required=True, help= "Path to positive dataset")
    parser.add_argument("-n", "--in_negative_data_filepath", required=True, help= "Path to negative dataset")
    parser.add_argument("-d", "--output_path", required=True, help= "Path to output directory")
    args = parser.parse_args()

    X, y = rl.read_pos_neg_data(args.in_positive_data_filepath, args.in_negative_data_filepath)
    # print(X.info())
    # print(y)
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=1)

    automl = autosklearn.classification.AutoSklearnClassifier(
                                                time_left_for_this_task=360, # increas means better mdel can be found
                                                per_run_time_limit=36, # default=1/10 of time_left_for_this_task
                                                initial_configurations_via_metalearning=25, # Disable if the hyperparameter optimization algorithm should start from scratch.
                                                ensemble_size =50, # Number of models from libraries of models
                                                ensemble_nbest=1, # only consider n best models to build ensemble
                                                max_models_on_disc=50, # n models keept on disc
                                                seed=1, # used to seed SMAC (outpufile name)
                                                memory_limit=3072, # MB memory limit for fitting the algos
                                                resampling_strategy='cv', # not sure how to set default holdout’ (train_size default=0.67)
                                                resampling_strategy_arguments={'folds': 10},
                                                metric=autosklearn.metrics.f1,
                                                tmp_folder=args.output_path + '/autosklearn_classification_example_tmp')
    automl.fit(X_train, y_train)
    # dataset_name='RRI_data'-> only for meta data


    df = automl.leaderboard(detailed=True)
    # need if cv
    # return a representation of the final ensemble found by auto-sklearn.
    print("Final ensamble\n:",automl.show_models())
    print("Final params\n:",automl.get_params())
    print("Stasitstic\n:",automl.sprint_statistics())

    print('Before re-fit')
    predictions = automl.predict(X_test)
    print("Accuracy score CV", sklearn.metrics.accuracy_score(y_test, predictions))


    print('After re-fit')
    automl.refit(X_train.copy(), y_train.copy())
    predictions = automl.predict(X_test)
    print("Accuracy score CV", sklearn.metrics.accuracy_score(y_test, predictions))

    # Return the underlying estimator object
    #model = automl.get_estimator()



    model_path = args.output_path + "/model.obj"
    print('output model: ', model_path)
    model_handle = open(model_path,"wb")
    pickle.dump(automl,model_handle)
    model_handle.close()



if __name__ == '__main__':
    main()
