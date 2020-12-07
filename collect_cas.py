from Bio import SeqIO, SearchIO
import json
import argparse
import os

# Processes putative Cas sequences from CRISPRCasFinder into an improved summary file indexed by contig.

# Parse command line arguments
parser = argparse.ArgumentParser(description='Processes putative Cas sequences from CRISPRCasFinder into an improved summary file.')
parser.add_argument('--result', action='store', dest='json_file_path', help='path to folder containing results.json from CCF')
parser.add_argument('--out', action='store', dest='out_path', help='where to store Cas output file, Cas_summary.tsv')

args = parser.parse_args()

# Parse CRISPRCasFinder results into a Python dictionary
json_dict = []
with open(args.json_file_path, 'r') as file:
    json_dict = json.load(file)
    file.close()
            
# Only useful data is contained in sequences
json_seq_list = json_dict['Sequences']

# Write a row for each cas gene detected, associated with a type and contig_id
with open(os.path.join(args.out_path, 'Cas_summary.tsv'), 'w') as f:
    for seq in json_seq_list:
        id_ = seq["Id"]
        for cas_cluster in seq['Cas']:
            cas_type = cas_cluster['Type']
            for gene in cas_cluster['Genes']:
                gene_subtype = gene['Sub_type']
                f.write(f"{id_}\t{cas_type}\t{gene_subtype}\n")
