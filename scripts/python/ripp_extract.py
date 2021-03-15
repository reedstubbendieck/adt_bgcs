#!/usr/bin/python3

from nlpprecursor.classification.data import DatasetGenerator as CDG
from nlpprecursor.annotation.data import DatasetGenerator as ADG
from pathlib import Path
import argparse, csv, json, nlpprecursor, os, sys
import pandas as pd

parser = argparse.ArgumentParser(
        description='''Use NLPPrecursor from the Magarvey lab to identify RiPP precursor peptides from a fasta file''')

# add the parser arguments
parser.add_argument("input", type=str, help="path/to/input/fasta.faa")
parser.add_argument("output", type=str, help="path/to/output/dir/")
parser.add_argument("models", type=str, help="path/to/models/")

# This allows for backwards compatibility of the pickled models.
sys.modules["protai"] = nlpprecursor

# parse positional arguments
args=parser.parse_args()

models_dir = Path(args.models)
input_path = args.input
output_path = args.output

# set models
class_model_dir = models_dir / "classification"
class_model_path = class_model_dir / "model.p"
class_vocab_path = class_model_dir / "vocab.pkl"
annot_model_dir = models_dir / "annotation"
annot_model_path = annot_model_dir / "model.p"
annot_vocab_path = annot_model_dir / "vocab.pkl"

# parse the input fasta file
fasta_dict = {}

with open(input_path) as n:
    for line in n:
        if line.startswith(">"):
            sequence_name_pre = line.rstrip().lstrip(">")
            sequence_name = sequence_name_pre.split(" ")[0]
        else:
            seq = fasta_dict.setdefault(sequence_name, "")
            fasta_dict[sequence_name] = seq + line.rstrip()

sequences = []

for sequence in fasta_dict:
    seq_dict = {"sequence":fasta_dict[sequence], "name":sequence}
    sequences.append(seq_dict)

class_predictions = CDG.predict(class_model_path, class_vocab_path, sequences)
cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, sequences)

ripp_class_dict = json.loads(json.dumps(class_predictions, indent=4))
cleavage_predictions_dict = json.loads(json.dumps(cleavage_predictions, indent=4))

# open the output_class text file
output_class = open(output_path+"ripp_classes.tsv", "w+")

## write the header line
output_class.write("orf\tclass\tscore\n")

## parse the ripp_class_dict to pull out the class predictions and scores
for i in range(0,len(ripp_class_dict)):
    name = ripp_class_dict[i]['name']
    for j in range(0,len(ripp_class_dict[i]['class_predictions'])):
        pred_dict = {}
        pred_dict[ripp_class_dict[i]['class_predictions'][j]['class']] = ripp_class_dict[i]['class_predictions'][j]['score']
        for class_type in pred_dict:
            output_class.write((name+"\t"+class_type+"\t"+str(pred_dict[class_type]))+"\n")

output_class.close()

# open the output_ripp_seq file
output_ripp_seq = open(output_path+"ripp_seq.tsv", "w+")

## write the header line
output_ripp_seq.write("orf\tsequence\tstatus\n")

## parse the cleavage_predictions_dict to pull out the core peptide sequences and statistics
for i in range(0,len(cleavage_predictions_dict)):
    name = cleavage_predictions_dict[i]['name']
    seq = cleavage_predictions_dict[i]['cleavage_prediction']['sequence']
    success = cleavage_predictions_dict[i]['cleavage_prediction']['status']
    output_ripp_seq.write(name+"\t"+seq+"\t"+success+"\n")

output_ripp_seq.close()

# join the two sequence tables and filter out low quality peptides

## load in the two output files from above into pandas dfs and merges them into a single dataframe
class_frame = pd.read_csv(output_path+"ripp_classes.tsv", sep="\t")
ripp_frame = pd.read_csv(output_path+"ripp_seq.tsv", sep="\t")
merged_frame = pd.merge(class_frame,ripp_frame,on='orf',how='left')

### filter out orfs that are of class NONRIPP, fail to yield a predicted structure, and possess a low score
merged_frame = merged_frame[(merged_frame['class'] != 'NONRIPP') & (merged_frame['status'] == 'success') & (merged_frame['score'] >= 0.75)]

### save the merged dataframe to file for later use as metadata
merged_frame.to_csv(output_path+"merged_ripp_frame.tsv", sep = "\t", header=True, index = False, index_label=False)

### extract the core peptide sequences into a dictionary then write to file
peptide_dict = dict(zip(merged_frame['orf'], merged_frame['sequence']))

# write the extracted sequences out to a file

output_fasta = open(output_path+"ripp_core_peptides.fasta", "w+")

for seq in peptide_dict:
    header = ">" + seq + "\n"
    output_fasta.write(header)
    output_fasta.write(peptide_dict[seq])
    output_fasta.write("\n")

output_fasta.close()