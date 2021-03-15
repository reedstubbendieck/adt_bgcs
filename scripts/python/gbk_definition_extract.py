#!/usr/bin/python3

# parse genbank files to extract definition lines and create table relating cluster to genomes

import argparse, csv, os
import pandas as pd
from pandas import Series

# parse command line arguments

parser = argparse.ArgumentParser(
        description='''Extract cluster number and SEQFID from .gbk files''')

# add the parser arguments
parser.add_argument("input_dir", type=str, help="path/to/dir/with/files/*.gbk")
parser.add_argument("meta_in", type=str, help="path/to/SEQFID_info.txt")
parser.add_argument("output", type=str, help="path/to/output/")

args=parser.parse_args()

input_path = args.input_dir
meta_in_path = args.meta_in
output_path = args.output

# extract cluster number and SEQFID from .gbk files

cluster_dict = {}
type_dict = {}

for file in os.listdir(input_path):
    if file.endswith(".gbk"):
        filename = file.split(".")[0]
        with open(input_path+file, 'r') as f:
            for line in f:
                if line.startswith("DEFINITION"):
                    organism=((line.split("DEFINITION  ")[1]).split("_")[0])
                    cluster_dict[filename] = organism
                if line.startswith("                     /product="):
                    type_dict[filename]=(line.split("\"")[1])

# associate metadata with clusters
## import the eHOMD SEQFID_info.txt file and extract the genus, species, and strain information for each cluster

metadata = pd.read_csv(meta_in_path, sep="\t")

genus_dict = {}
species_dict = {}
strain_dict = {}

for cluster in cluster_dict:
    x = (metadata[metadata['SEQFID']==cluster_dict[cluster]]).iloc[0]
    genus_dict[cluster]=x.iloc[4]
    species_dict[cluster]=x.iloc[5]
    strain_dict[cluster]=str(x.iloc[6])


# associate clusters with BGC type

with open(output_path, 'w', newline='') as f_output:
    f_output.write("cluster\tSEQFID\tgenus\tspecies\tstrain\n")
    for item in cluster_dict:
        x = item+"\t"+cluster_dict[item]+"\t"+genus_dict[item]+"\t"+species_dict[item]+"\t"+strain_dict[item]+"\n"
        f_output.write(x)
