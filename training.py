#!/usr/bin/env python
import pandas as pd
import sys
import argparse
import rrieval.lib as rl

def main():
    parser = argparse.ArgumentParser(description='Trains models for ')
    parser.add_argument("-i", "--in_path", action="store", dest="input_path",
                        required=True,
                        help= "path to folder storing all input data")
    parser.add_argument("-d", "--output_path", action="store", dest="output_path",
                        required=True,
                        help= "path output reposetory")

    args = parser.parse_args()

if __name__ == '__main__':
    main()
