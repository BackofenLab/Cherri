import sys
import argparse
import rrieval.lib as rl

def main():
    parser = argparse.ArgumentParser(description='Trains models for RRIeval')
    parser.add_argument("-ip", "--in_positive_data_filepath", required=True, help= "Path to positive dataset")
    parser.add_argument("-in", "--in_negative_data_filepath", required=True, help= "Path to negative dataset")
    parser.add_argument("-d", "--output_path", required=True, help= "Path to output directory")
    parser.add_argument("-n", "--name", required=True, help= "Name of base model")
    args = parser.parse_args()

    # compute feature importance!

    rl.base_model(args.in_positive_data_filepath,args.in_negative_data_filepath,args.output_path,args.name)

if __name__ == '__main__':
    main()
