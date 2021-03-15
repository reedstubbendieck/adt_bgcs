#!/usr/bin/env python
# coding: utf-8

# extract top hits from HMM scan output and filter out non-biosynthethic ORFs commonly found in BGCs

import argparse, csv, os, subprocess
import pandas as pd

parser = argparse.ArgumentParser(
        description='''Processes .gbk files from antiSMASH and generates nucleotide sequences without non-BGC orfs''')

# add the parser arguments
parser.add_argument("input_dir", type=str, help="path/to/dir/with/inputs/")
parser.add_argument("exclude_list", type=str, help="path/to/exclude/list.txt")
parser.add_argument("hmm", type=str, help="path/to/Pfam-A.hmm")

args=parser.parse_args()

input_path = args.input_dir
exclude_path = args.exclude_list
hmm_path = args.hmm

excludeList = [line.rstrip('\n') for line in open(exclude_path)]

# function for processing .gbk files into .faa and .fnas using genbank_to_fasta.py

def asProcessor(input_gbk):
    prefix = os.path.split(input_gbk)[1].split(".gbk")[0]
    if not os.path.isfile(input_path+prefix+".fna"):
        fasta_convert_nt_call = "genbank_to_fasta.py -i "+input_gbk+" --out_file "+prefix+".fna -m genbank -s nt -f CDS -q locus_tag"
        subprocess.call(fasta_convert_nt_call, shell=True)
    else:
        print(prefix+".fna exists")
    if not os.path.isfile(input_path+prefix+".faa"):
        fasta_convert_aa_call = "genbank_to_fasta.py -i "+input_gbk+" --out_file "+prefix+".faa -m genbank -s aa -f CDS -q locus_tag"
        subprocess.call(fasta_convert_aa_call, shell=True)
    else:
        print(prefix+".faa exists")

# function for running hmmscan on .faa files

def hmmscan(input_faa):
    prefix = os.path.split(input_faa)[1].split(".faa")[0]
    if not os.path.isfile(input_path+prefix+"_hmm"):
        hmmscan_call = "hmmscan --tblout "+input_path+prefix+"_hmm --acc "+hmm_path+" "+input_faa
        subprocess.call(hmmscan_call, shell=True)
    else:
        print(prefix+"_hmm exists")

# function for reading hmmscan output into df and processing

def hmm_parser(input_hmm):
    # parse the hmmscan output and pull out relevant columns
    global pass_check
    pass_check = 1
    try:
        hmm_out = pd.read_csv(input_hmm, comment='#', delimiter=r"\s+", header=None, usecols=[0,1,2,3,4,5])
        # remove the decimals from PFAM domains
        hmm_out[1] = hmm_out[1].str.replace('\..*', '')
        # group the hmm output by query and find best hit
        hmm_out.groupby([2], sort=False)[4].min()
        idx = hmm_out.groupby([2])[4].transform(min) == hmm_out[4]
        hmm_out = hmm_out[idx]
        # add the excluded hits into counter dictionary
        global excluded_orf_count
        hmm_filtered_out = hmm_out[hmm_out[1].isin(excludeList)]
        excluded_orfs = hmm_filtered_out[1].tolist()
        for orf in excluded_orfs:
            if orf in excluded_orf_count:
                excluded_orf_count[orf] = excluded_orf_count[orf] + 1
            else:
                excluded_orf_count[orf] = 1
        # only return hits that are not in the excluded list
        hmm_filtered_in = hmm_out[~hmm_out[1].isin(excludeList)]
        included_orfs = hmm_filtered_in[2].tolist()
        return(included_orfs)
    # raise exception for hmm files that contain no predicted domains
    except:
        print(input_hmm+" has an error, please check manually")
        exception_list.append(input_hmm)
        pass_check = 0
        pass

# function for generating files that just include the biosynthethic ORFs

def fasta_handler(input_fna,x):
    fasta_dict = {}
    with open(input_fna) as n:
        for line in n:
            if line.startswith(">"):
                sequence_name_pre = line.rstrip().lstrip(">")
                sequence_name = sequence_name_pre.split(" ")[0]
            else:
                seq = fasta_dict.setdefault(sequence_name, "")
                fasta_dict[sequence_name] = seq + line.rstrip()
    global in_prefix
    if not os.path.isfile(input_path+in_prefix+"_biosynthetic.fna"):
        # write file with only biosynthetic ORFs
        with open (input_path+in_prefix+"_biosynthetic.fna", "w") as u:
            for sequence in fasta_dict:
                if sequence in x:
                    header = ">" + in_prefix + "_" + sequence + "\n"
                    u.write(header)
                    u.write(fasta_dict[sequence])
                    u.write("\n")
    else:
        print(in_prefix+"_biosynthetic.fna exists")

# initialize the excluded orf count dictionary
excluded_orf_count = {}

# initialize the exception list

exception_list = []

for file in os.listdir(input_path):
        if file.endswith(".gbk"):
            in_prefix = file.split(".gbk")[0]
            print("Processing "+in_prefix)
            asProcessor(input_path+file)
            hmmscan(input_path+in_prefix+".faa")
            x=hmm_parser(input_path+in_prefix+"_hmm")
            if pass_check == 1:
                fasta_handler(input_path+in_prefix+".fna",x)

# write the global excluded orf count to a .csv file
(pd.DataFrame.from_dict(data=excluded_orf_count, orient='index')
   .to_csv('./derivedData/excluded_orf_count.csv', header=False))

# check these manually

print("The following hmm files were invalid. Please check manually. Will move the associated files to the bad_clusters/ directory")
for p in exception_list:
    print(p)

if not os.path.exists(input_path+"bad_clusters/"):
    os.makedirs(input_path+"bad_clusters/")
    
for cluster in exception_list:
    move_prefix = (os.path.basename(cluster).split("_hmm")[0])
    for file in os.listdir(input_path):
        if file.startswith(move_prefix):
            os.rename(input_path+file, input_path+"bad_clusters/"+file)