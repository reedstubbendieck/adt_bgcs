#!/usr/bin/env python
# processes reads with fastp then aligns to index
# modified to deal with files from hmp stored as .tar.bz2

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
    if file.endswith(".tar.bz2"):
        # identify unique prefix for directory naming
        file_prefix=(hmp_cart[hmp_cart['urls'].str.contains(file)]['file_id'].tolist()[0])
        print(file_prefix)
        if not os.path.exists(output_path+file_prefix):
            os.makedirs(output_path+file_prefix)
        # check if kallisto has already been run, if so moves on and deltes the .tar.bz2 file
        if os.path.isfile(output_path+file_prefix+"/abundance.tsv"):
            print(output_path+file_prefix+"/abundance.tsv exists")
            del_call_1 = "rm -rf "+input_path+file
            subprocess.call(del_call_1, shell=True)
        else:
            # uncompress the reads
            tar_call = "pbzip2 -dcvm2000 "+input_path+file+" | tar -C "+input_path+" -vx"
            subprocess.call(tar_call, shell=True)
            # need to move files from directory in *.tar.bz2 file to current directory
            mv_call = 'find . -name "*1.fastq" -exec mv "{}" '+input_path+' \;'
            mv_call2 = 'find . -name "*2.fastq" -exec mv "{}" '+input_path+' \;'
            # use the del_call to delete the directories that are generated during the untar
            del_call = "rm -rf "+input_path+file.split(".tar.bz2")[0]
            del_call_bz = "rm -f "+input_path+file
            subprocess.call(mv_call, shell=True)
            subprocess.call(mv_call2, shell=True)
            subprocess.call(del_call, shell=True)
            subprocess.call(del_call_bz, shell=True)
            for filegz in os.listdir(input_path):
                # need to run fastp on all of the files from the tar archives. The following uncompresses the tar files and then runs fastp in PE mode on corresponding reads
                if filegz.endswith("1.fastq"):
                    read_prefix = filegz.split(".1.fastq")[0]
                    fastp_call = "fastp -i "+input_path+filegz+" -I "+input_path+read_prefix+".2.fastq -o "+input_path+read_prefix+"_R1.out.fastq.gz -O "+input_path+read_prefix+"_R2.out.fastq.gz -h "+input_path+read_prefix+".html -j "+input_path+read_prefix+".json -w "+threads_num
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
                            # clean up the uncompressed reads and tar.bz2 files
                            del_call_3 = "rm -rf "+input_path+file.split(".tar.bz2")[0]+"*"
                            subprocess.call(del_call_3, shell=True)
