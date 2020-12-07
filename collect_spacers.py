from Bio import SeqIO, SearchIO
import json
import argparse
import os

# Processes spacer sequences from CRISPRCasFinder into an indexed form which can be BLASTed.

# Parse command line arguments
parser = argparse.ArgumentParser(description='Filter contigs for those containing CRISPR/Cas systems')
parser.add_argument('--result', action='store', dest='json_file_path', help='path to folder containing results.json from CCF')
parser.add_argument('--out', action='store', dest='out_path', help='where to store spacer output file, labelled_spacers.fa')
args = parser.parse_args()

# Parse CRISPRCasFinder results into a Python dictionary
json_dict = []
with open(args.json_file_path, 'r') as file:
    json_dict = json.load(file)
    file.close()
            
# Only useful data is contained in sequences
json_seq_list = json_dict['Sequences']

# Write a unique name and spacer sequence to a new file
with open(os.path.join(args.out_path, 'labelled_spacers.fa'), 'w') as f:
    for seq in json_seq_list:
        id_ = seq["Id"]
        crispr_no = 0
        spacer_no = 0
        for crispr in seq['Crisprs']:
            regions = crispr['Regions']
            for region in regions:
                if region['Type'] == "Spacer":
                    f.write(">"+"_".join([str(id_), str(crispr_no), str(spacer_no)]))
                    f.write('\n'+region['Sequence']+'\n')
                    spacer_no += 1
            crispr_no += 1
