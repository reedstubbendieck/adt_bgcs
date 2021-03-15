#!/usr/bin/python3

# rename massive fasta file headers to be antiSMASH compliant

import argparse

# parse command line arguments

parser = argparse.ArgumentParser(
        description='''trim contig names from eHOMD genome assemblies''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/genomes.fna")
parser.add_argument("output", type=str, help="path/to/output/genomes.fna/")

args=parser.parse_args()

input_path = args.input
output_path = args.output

fasta_dict = {}
sequence_header_list = []
sequence_header_count = {}

with open(input_path) as n:
    for line in n:
        if line.startswith(">"):
            sequence_name = line.rstrip().lstrip(">")
            header = sequence_name.split("_")[0]
            if header not in sequence_header_list:
                sequence_header_list.append(header)
                sequence_header_count[header] = 1
            else:
                sequence_header_count[header] = sequence_header_count[header] + 1
            sequence_head = header+"_c"+str(sequence_header_count[header])
        else:
            seq = fasta_dict.setdefault(sequence_head, "")
            fasta_dict[sequence_head] = seq + line.rstrip()

# write the codon alignments out to a file

output_file = open(output_path, "w+")

for sequence in fasta_dict:
    header = ">" + sequence + "\n"
    output_file.write(header)
    output_file.write(fasta_dict[sequence])
    output_file.write("\n")

output_file.close()
