#!/usr/bin/env python
# kallisto_quant_processing.py

import argparse, os, subprocess
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(
        description='''combines the metagenome read abundance.tsv files from kallisto into a matrix''')

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

tsv_list = []
data_frames = []

for entry in Path(input_path).iterdir():
    # check if entry in path is a directory
    if entry.is_dir():
        # checks if the abundance.tsv file exists in directory
        if os.path.isfile(input_path+entry.name+"/abundance.tsv"):
            tsv_list.append(input_path+entry.name+"/abundance.tsv")

for tsv in tsv_list:
    frame = pd.read_csv(tsv, delimiter="\t")
    frame['file_id'] = os.path.dirname(tsv).rsplit('/', 1)[1]
    data_frames.append(frame)
 

bigframe = pd.concat(data_frames, ignore_index=True)

# generates a cluster summary from the bigframe

## drop unncessary columns from the data frame
bigframe = bigframe.drop('length', axis=1)
bigframe = bigframe.drop('eff_length', axis=1)
bigframe = bigframe.drop('tpm', axis=1)

## remove the SEQF1003_c2. header from the target_id column
bigframe['target_id'] = bigframe['target_id'].replace({'SEQF1003_c2.':''}, regex=True)

bigframe[['cluster','orf']] = bigframe['target_id'].str.split("_", 1, expand=True) 

## clean up the data frame
bigframe = bigframe.drop('target_id', axis=1)
bigframe = bigframe.drop('orf', axis=1)

## summarize by file_id and cluster
bigframe = bigframe.groupby(['file_id','cluster']).agg({'est_counts': 'sum'}).reset_index()

## sava the data frame
bigframe.to_csv(output_path+"hmp_adt_kallisto_counts.tsv", sep = "\t", header=True, index = False, index_label=False)