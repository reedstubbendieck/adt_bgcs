#!/usr/bin/env python

import argparse, csv, os, json
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(
        description='''extract read counts''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/dir/")
parser.add_argument("output", type=str, help="path/to/output/dir/")

args=parser.parse_args()

input_path = args.input
output_path = args.output

# make the output directory if it doesn't exist

if not os.path.exists(output_path):
    os.makedirs(output_path)

# generate the list of abundance.tsv files

read_num_dict = {}

for entry in Path(input_path).iterdir():
    # check if entry in path is a directory
    if entry.is_dir():
        # checks if the abundance.tsv file exists in directory
        if os.path.isfile(input_path+entry.name+"/run_info.json"):
            f = open(input_path+entry.name+"/run_info.json")
            json_file = json.load(f)
            read_num_dict[entry.name] = (json_file['n_processed'])
            
# write the dict to a .csv file
(pd.DataFrame.from_dict(data=read_num_dict, orient='index')
   .to_csv(output_path+'read_count.csv', header=False))