#!/usr/bin/env python
import sys
import argparse
import rrieval.lib as rl

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-i", "--in_positive_data_filepath", required=True, help= "Path to positive dataset")
    parser.add_argument("-n", "--in_negative_data_filepath", required=True, help= "Path to negative dataset")
    parser.add_argument("-d", "--output_path", required=True, help= "Path to output directory")
    args = parser.parse_args()

    # compute feature importance!



    rl.param_optimize(args.in_positive_data_filepath,args.in_positive_data_filepath,args.output_path)

if __name__ == '__main__':
    main()
