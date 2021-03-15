#!/usr/bin/env python
# metaphlan_table_processing.py

import argparse, os, subprocess
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(
        description='''combines the individual metaphlan output tsvs into a single file for downstream processing''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/dir/")
parser.add_argument("output", type=str, help="path/to/output/dir/")
parser.add_argument("level", type=str, help="taxa level: genus/species")

args=parser.parse_args()

input_path = args.input
output_path = args.output
taxa_level = args.level

# make the output directory if it doesn't exist

if not os.path.exists(output_path):
    os.makedirs(output_path)
    
# generate the list of _classified files

tsv_list = []
data_frames = []

for file in os.listdir(input_path):
    if file.endswith("_classified"):
        tsv_list.append(os.path.join(input_path, file))

for tsv in tsv_list:
    frame = pd.read_csv(tsv, delimiter="\t", skiprows=4, header = None)
    file_id = (tsv.rsplit('/', 1)[1]).split("_")[0]
    frame.columns = [taxa_level, "NCBI_Tax_ID", "rel_abundance", "drop"]
    # remove the NCBI_Tax_ID and drop columns
    frame.drop("drop", axis=1, inplace=True)
    frame.drop("NCBI_Tax_ID", axis=1, inplace=True)
    # add the file_id column
    frame['file_id'] = file_id
    # remove the prefix from each taxa
    if taxa_level == "genus":
        frame[taxa_level] = frame[taxa_level].str.replace("g__", "")
    elif taxa_level == "species":
        frame[taxa_level] = frame[taxa_level].str.replace("s__", "")
    data_frames.append(frame)
    
bigframe = pd.concat(data_frames, ignore_index=True)

bigframe.to_csv(output_path+"hmp_adt_metaphlan_"+taxa_level+".tsv", sep = "\t", header=True, index = False, index_label=False)