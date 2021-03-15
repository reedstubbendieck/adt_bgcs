#!/usr/bin/env python
# processes reads with fastp then classifies with metaphlan

import argparse, os, subprocess
import pandas as pd

parser = argparse.ArgumentParser(
        description='''processes reads with fastp and classifies using metaphlan''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/dir/")
parser.add_argument("manifest", type=str, help="path/to/manifest/")
parser.add_argument("output", type=str, help="path/to/output/dir/")
parser.add_argument("level", type=str, help="taxa level: (g)enus or (s)pecies")
parser.add_argument("threads", type=str, help="number of threads")

args=parser.parse_args()

input_path = args.input
manifest_path = args.manifest
output_path = args.output
taxa_level = args.level
threads_num = args.threads

# going to use the file_id in the hmp cart to name the outputs

hmp_cart = pd.read_csv(manifest_path, delimiter="\t")

# uncompresses the raw reads from the iHMP
for file in os.listdir(input_path):
    if file.endswith(".fastq.tar"):
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
                fastp_call = "fastp -i "+input_path+filegz+" -I "+input_path+read_prefix+"_R2.fastq.gz -o "+output_path+file_prefix+"/"+read_prefix+"_R1.out.fastq.gz -O "+output_path+file_prefix+"/"+read_prefix+"_R2.out.fastq.gz -h "+output_path+file_prefix+"/"+read_prefix+".html -j "+output_path+file_prefix+"/"+read_prefix+".json -w "+threads_num
                subprocess.call(fastp_call, shell=True)
        # the following block finds all of the R1.out.fastq.gz files and then puts their corresponding PE in a list to feed into kallisto in the right order
        files = []
        for fastpfilegz in os.listdir(output_path+file_prefix+"/"):
            if fastpfilegz.endswith("R1.out.fastq.gz"):
                files.append(output_path+file_prefix+"/"+fastpfilegz)
                paired_file = output_path+file_prefix+"/"+fastpfilegz.split("R1.out.fastq.gz")[0]+"R2.out.fastq.gz"
                files.append(paired_file)
        # need to extract the forward and reverse reads into separate lists
        forward_files = files[0::2]
        reverse_files = files[1::2]
        forward_input = ' '.join([str(v) for v in forward_files])
        reverse_input = ' '.join([str(v) for v in reverse_files])
        # cat the files together
        forward_cat = "cat "+forward_input+" > "+output_path+file_prefix+"/"+"forward.fastq.gz"
        reverse_cat = "cat "+reverse_input+" > "+output_path+file_prefix+"/"+"reverse.fastq.gz"
        subprocess.call(forward_cat, shell=True)
        subprocess.call(reverse_cat, shell=True)
        # calls metaphlan
        metaphlan_call = "metaphlan "+output_path+file_prefix+"/"+"forward.fastq.gz,"+output_path+file_prefix+"/"+"reverse.fastq.gz --nproc "+threads_num+" --ignore_eukaryotes --ignore_archaea --tax_lev "+taxa_level+" --force --bowtie2out "+output_path+file_prefix+"/"+read_prefix+".bowtie --input_type fastq > "+output_path+file_prefix+"_classified"
        subprocess.call(metaphlan_call, shell=True)
        # clean up the uncompressed read
        cleanup_call = "rm -f "+input_path+file.split(".raw.fastq.tar")[0]+"*"
        cleanup_call2 = "rm -rf "+output_path+file_prefix+"/"
        subprocess.call(cleanup_call, shell=True)
        subprocess.call(cleanup_call2, shell=True)