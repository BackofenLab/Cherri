#!/usr/bin/env python
import sys
import argparse
import rrieval.lib as rl
import sklearn.datasets
import sklearn.metrics

import autosklearn.classification

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
                                                time_left_for_this_task=2000,
                                                per_run_time_limit=300,
                                                tmp_folder=args.output_path + '/autosklearn_classification_example_tmp',
                                                output_folder=args.output_path + '/autosklearn_classification_example_out')
    automl.fit(X_train, y_train, dataset_name='RRI_data')
    # return a representation of the final ensemble found by auto-sklearn.
    print(automl.show_models())

    predictions = automl.predict(X_test)
    print("Accuracy score:", sklearn.metrics.accuracy_score(y_test, predictions))

    print(automl.sprint_statistics())



if __name__ == '__main__':
    main()
