#!/usr/bin/python3

# parse genbank files to extract definition lines and create table relating cluster to genomes

import argparse, csv, os
import pandas as pd
from collections import defaultdict

# parse command line arguments

parser = argparse.ArgumentParser(
        description='''Extract t1pks and nrps domain information from .gbk files''')

# add the parser arguments
parser.add_argument("input_dir", type=str, help="path/to/dir/with/files/*.gbk")
parser.add_argument("output", type=str, help="path/to/output.tsv")

args=parser.parse_args()

input_path = args.input_dir
output_path = args.output

type_dict = {}
size_dict = {}
process_list = []
edge_dict = {}
domain_dict = defaultdict(dict) # initialize dictionary to store domain counts for each cluster

# extracts the product type from .gbk files
for file in os.listdir(input_path):
    if file.endswith(".gbk"):
        filename = file.split(".gbk")[0]
        with open(input_path+file, 'r') as f:
            for line in f:
                if line.startswith("                     /product="):
                    type_dict[filename]=(line.split("\"")[1])
                    process_list.append(file)
                if line.startswith("     cluster         1.."):
                    size_dict[filename]=((line.split("..")[1]).rstrip())
                    
# processes clusters to extract contig_edge domain counts
for file in process_list:
    with open(input_path+file, 'r') as f:
            filename=file.split(".gbk")[0]
            domain_list = [] # list of domains for each .gbk
            for line in f:
                if line.startswith("                     /contig_edge="):
                    edge_dict[filename] = (line.split("\"")[1])
                if line.startswith("                     /domain="):
                    domain_type=(line.split("\"")[1])
                    # check if this domain type has already been counted
                    if domain_type in domain_list:
                        domain_dict[filename][domain_type] = domain_dict[filename][domain_type] + 1
                    else:
                        # append the domain type to the list and initialize the counter
                        domain_list.append(domain_type)
                        domain_dict[filename][domain_type] = 1

# save the data to output files
## first file is the cluster_info
with open(output_path+"eHOMD_cluster_info.tsv", 'w', newline='') as f_output:
    f_output.write("cluster\ttype\tsize\tcontig_edge\n")
    for cluster in process_list:
        cluster_name = cluster.split(".gbk")[0]
        x = cluster_name+"\t"+type_dict[cluster_name]+"\t"+size_dict[cluster_name]+"\t"+edge_dict[cluster_name]+"\n"
        f_output.write(x)

## second file is domain counts per cluster
with open(output_path, 'w', newline='') as f_output:
    f_output.write("cluster\tdomain\tcounts\n")
    for cluster in domain_dict:
        for domain_type in domain_dict[cluster]:
            x = cluster+"\t"+domain_type+"\t"+str(domain_dict[cluster][domain_type])+"\n"
            f_output.write(x)
    