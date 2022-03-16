#!/usr/bin/env python
import sys
import argparse
import rrieval.lib as rl

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--in_data_filepath", required=True, help= "Path to input dataset")
    parser.add_argument("-m", "--in_model_filepath", required=True, help= "Path to trained model")
    parser.add_argument("-d", "--output_path", required=True, help= "Path to output directory")
    args = parser.parse_args()

    y_pred =  rl.classify(args.in_data_filepath,args.in_model_filepath,args.output_path)

    pos = 0
    neg = 0

    for i in y_pred:
        if i == 1:
            pos +=1
        elif i == 0:
            neg +=1
    print('pos: %i \nneg: %i'%(pos, neg))

if __name__ == '__main__':
    main()
