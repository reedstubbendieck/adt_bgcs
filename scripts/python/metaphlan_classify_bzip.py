#!/usr/bin/env python
# processes reads with fastp then uses metaphlan to classify them

import argparse, os, subprocess, tarfile
import pandas as pd

parser = argparse.ArgumentParser(
        description='''processes reads with fastp and classifies reads with metaphlan''')

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

# make the output directory if it doesn't exist

if not os.path.exists(output_path):
    os.makedirs(output_path)

# going to use the file_id in the hmp cart to name the outputs

hmp_cart = pd.read_csv(manifest_path, delimiter="\t")

# the following code splits the manifest file into individual manifest files with one row per file named after the file_id

for i, g in hmp_cart.groupby("file_id"):
    g.to_csv(output_path+'{}.hmp.manifest'.format(i), sep = "\t", header=True, index = False, index_label=False)
    
# generate list of file_id to process

process_list = hmp_cart['file_id'].tolist()

# loop through the manifest files to download from iHMP using portal_client

for file in os.listdir(output_path):
    if file.endswith("hmp.manifest"):
        hmp_file_id = file.split(".hmp.manifest")[0]
        # check to see if the file has been run previously
        if os.path.isfile(output_path+hmp_file_id+"_classified"):
            print(output_path+hmp_file_id+"_classified exists, skipping and deleting manifest file")
            manifest_del_call = "rm -rf "+output_path+file
            subprocess.call(manifest_del_call, shell=True)
        else:
            if hmp_file_id in process_list:
                # uncompress the downloaded reads, uses the urls column to specify the file name of the compressed reads
                compressed_reads = ((hmp_cart[hmp_cart['file_id'].str.contains(hmp_file_id)]['urls'].tolist()[0]).split(",")[0]).rsplit('/', 1)[1]
                # check to make sure file actually downloaded and creates directories to store the metaphlan output and processed reads
                if not os.path.exists(output_path+"reads/"+hmp_file_id):
                    os.makedirs(output_path+"reads/"+hmp_file_id)
                # extract files from the downloaded compressed reads
                compressed_path = input_path+compressed_reads
                if compressed_path.endswith(".tar"):
                    tar = tarfile.open(input_path+compressed_reads)
                    tar.extractall(output_path+"reads/"+hmp_file_id+"/")
                else:
                    tar = tarfile.open(input_path+compressed_reads, "r:*")
                    for item in tar:
                        if item.name.endswith(".1.fastq") or item.name.endswith(".2.fastq"):
                            item.name = os.path.basename(item.name)
                            tar.extract(item, output_path+"reads/"+hmp_file_id+"/")
                tar.close()
                for filegz in os.listdir(output_path+"reads/"+hmp_file_id+"/"):
                    # need to run fastp on all of the files from the tar archives. The following runs fastp in PE mode on corresponding reads
                    if filegz.endswith("1.fastq"):
                        read_prefix = filegz.split(".1.fastq")[0]
                        fastp_call = "fastp -i "+output_path+"reads/"+hmp_file_id+"/"+filegz+" -I "+output_path+"reads/"+hmp_file_id+"/"+read_prefix+".2.fastq -o "+output_path+"reads/"+hmp_file_id+"/"+read_prefix+"_R1.out.fastq.gz -O "+output_path+"reads/"+hmp_file_id+"/"+read_prefix+"_R2.out.fastq.gz -h "+output_path+"reads/"+hmp_file_id+"/"+read_prefix+".html -j "+output_path+"reads/"+hmp_file_id+"/"+read_prefix+".json -w "+threads_num
                        subprocess.call(fastp_call, shell=True)
                        # the following block finds all of the R1.out.fastq.gz files, must combine all forward and reverse reads into a single file for each direction
                        files = []
                        for fastpfilegz in os.listdir(output_path+"reads/"+hmp_file_id+"/"):
                            if fastpfilegz.endswith("R1.out.fastq.gz"):
                                files.append(output_path+"reads/"+hmp_file_id+"/"+fastpfilegz)
                                paired_file = output_path+"reads/"+hmp_file_id+"/"+fastpfilegz.split("R1.out.fastq.gz")[0]+"R2.out.fastq.gz"
                                files.append(paired_file)
                                # need to extract the forward and reverse reads into separate lists
                                forward_files = files[0::2]
                                reverse_files = files[1::2]
                                forward_input = ' '.join([str(v) for v in forward_files])
                                reverse_input = ' '.join([str(v) for v in reverse_files])
                                # cat the files together
                                forward_cat = "cat "+forward_input+" > "+output_path+"reads/"+hmp_file_id+"/"+"forward.fastq.gz"
                                reverse_cat = "cat "+reverse_input+" > "+output_path+"reads/"+hmp_file_id+"/"+"reverse.fastq.gz"
                                subprocess.call(forward_cat, shell=True)
                                subprocess.call(reverse_cat, shell=True)
                                # calls metaphlan
                                metaphlan_call = "metaphlan "+output_path+"reads/"+hmp_file_id+"/"+"forward.fastq.gz,"+output_path+"reads/"+hmp_file_id+"/"+"reverse.fastq.gz --nproc "+threads_num+" --ignore_eukaryotes --ignore_archaea --tax_lev "+taxa_level+" --force --bowtie2out "+output_path+"reads/"+hmp_file_id+".bowtie --input_type fastq > "+output_path+hmp_file_id+"_classified"
                                subprocess.call(metaphlan_call, shell=True)
                                # clean up the uncompressed reads and tar.bz2 files
                                del_call = "rm -rf "+output_path+"reads/"+hmp_file_id+"/"
                                del_call2 = "rm -rf "+input_path+compressed_reads
                                del_call3 = "rm -rf "+output_path+file
                                subprocess.call(del_call, shell=True)
                                subprocess.call(del_call2, shell=True)
                                subprocess.call(del_call3, shell=True)