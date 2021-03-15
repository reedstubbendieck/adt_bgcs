#!/usr/bin/env python
# processes reads with fastp then aligns to index

import argparse, os, subprocess
import pandas as pd

parser = argparse.ArgumentParser(
        description='''processes reads with fastp and aligns to kallisto index''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/dir/")
parser.add_argument("manifest", type=str, help="path/to/manifest/")
parser.add_argument("index", type=str, help="path/to/kallisto/index")
parser.add_argument("output", type=str, help="path/to/output/dir/")
parser.add_argument("threads", type=str, help="number of threads")

args=parser.parse_args()

input_path = args.input
manifest_path = args.manifest
index_path = args.index
output_path = args.output
threads_num = args.threads

# going to use the file_id in the hmp cart to name the outputs

hmp_cart = pd.read_csv(manifest_path, delimiter="\t")

# uncompresses the raw reads from the iHMP
for file in os.listdir(input_path):
    if file.endswith(".fastq.tar"):
    # looks for matches in the hmp_cart urls column and returns fild_id for unique naming
        tar_call = "tar -C "+input_path+" -xvf "+input_path+file
        subprocess.call(tar_call, shell=True)
        # identify unique prefix for directory naming
        file_prefix=(hmp_cart[hmp_cart['urls'].str.contains(file.split(".tar")[0]+".tar")]['file_id'].tolist()[0])
        if not os.path.exists(output_path+file_prefix):
            os.makedirs(output_path+file_prefix)
        # need to run fastp on all of the files from the tar archives. The following uncompresses the tar files and then runs fastp in PE mode on corresponding reads
        for filegz in os.listdir(input_path):
            if filegz.endswith("_R1.fastq.gz"):
                read_prefix = filegz.split("_R1.fastq.gz")[0]
                fastp_call = "fastp -i "+input_path+filegz+" -I "+input_path+read_prefix+"_R2.fastq.gz -o "+input_path+read_prefix+"_R1.out.fastq.gz -O "+input_path+read_prefix+"_R2.out.fastq.gz -h "+input_path+read_prefix+".html -j "+input_path+read_prefix+".json -w "+threads_num
                subprocess.call(fastp_call, shell=True)
        # the following block finds all of the R1.out.fastq.gz files and then puts their corresponding PE in a list to feed into kallisto in the right order
        files = []
        for fastpfilegz in os.listdir(input_path):
            if fastpfilegz.endswith("R1.out.fastq.gz"):
                files.append(input_path+fastpfilegz)
                paired_file = input_path+fastpfilegz.split("R1.out.fastq.gz")[0]+"R2.out.fastq.gz"
                files.append(paired_file)
        kallisto_input = ' '.join([str(v) for v in files])
        # calls kallisto
        kallisto_call = "kallisto quant -i "+index_path+" -o "+output_path+file_prefix+"/ -t "+threads_num+" "+kallisto_input
        subprocess.call(kallisto_call, shell=True)
        # clean up the uncompressed read
        cleanup_call = "rm -f "+input_path+file.split(".raw.fastq.tar")[0]+"*"
        subprocess.call(cleanup_call, shell=True)
