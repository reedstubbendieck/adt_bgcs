#!/usr/bin/env python
# substitute the SEQFIDXXXX_c# in antiSMASH .gbk files with the cluster number before using BiG-SCAPE

import argparse, os, subprocess

parser = argparse.ArgumentParser(
        description='''substitute the SEQFIDXXXX_c# in antiSMASH .gbk files with the cluster number before using BiG-SCAPE''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/dir/")
parser.add_argument("output", type=str, help="path/to/output/dir/")

args=parser.parse_args()

input_path = args.input
output_path = args.output

# make the output directory if it doesn't exist

if not os.path.exists(output_path):
    os.makedirs(output_path)
	
for file in os.listdir(input_path):
    if file.endswith(".gbk"):
        cluster_prefix = file.split(".")[0]
        sed_call = "sed -r 's/SEQF[0-9]+_c[0-9]+/"+cluster_prefix+"/g' "+input_path+file+" > "+output_path+cluster_prefix+".gbk"
        subprocess.call(sed_call, shell=True)