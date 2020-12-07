from Bio import SeqIO, SearchIO
import argparse
import pandas as pd
import os

# Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table

parser = argparse.ArgumentParser(description='Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table')
parser.add_argument('--contigs', action='store', dest='contig_file_path', help='full path to original contigs fasta file')
parser.add_argument('--spacers', action='store', dest='spacer_file_path', help='full path to labelled_spacers.fa file')
parser.add_argument('--contigs_blast', action='store', dest='blast_result', help='full path to blast results tab from contig sequences')
parser.add_argument('--spacers_blast', action='store', dest='spacer_result', help='full path to blast results tab from spacer sequences')
#parser.add_argument('--result', action='store', dest='result_json', help='full path to CCF results.json file')
parser.add_argument('--out', action='store', dest='out_path', default='')

args=parser.parse_args()
if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)

# get all contig ids from contig file
contig_file = list(SeqIO.parse(args.contig_file_path, "fasta"))
all_contig_names = set([contig.id for contig in contig_file])

# parse contig blast results into a dataframe
df = pd.read_csv(args.blast_result, delimiter='\t', header=None) 
df.columns=['contig_id', 'accession', 'common_name', 'e_value']
# take only the top hit, by e-value. In the case of ties, this result is essentially random
df = df.drop_duplicates(subset='contig_id', keep='first')
# determine contigs with no blast hits
contigs_with_match = set(df['contig_id'])
unmatched_contigs = all_contig_names.difference(contigs_with_match)
# add unmatched contigs to dataframe, so that all contigs are represented
for unmatched in unmatched_contigs:
    df = df.append(pd.Series({'contig_id': unmatched, 
                              'accession': "no_match", 
                              'common_name': None,
                              'e_value': None}), ignore_index=True)

# get all spacer ids from spacer fle
spacer_file = list(SeqIO.parse(args.spacer_file_path, "fasta"))
all_spacer_names = set([spacer.id for spacer in spacer_file])

# parse blast results for spacers into dataframe
df_s = pd.read_csv(args.spacer_result, delimiter='\t', header=None)
df_s.columns=['contig_CRISPR_spacer', 'accession', 'common_name', 'e_value', 'bitscore', 'percent_ident']
# keep only the top hit for each spacer, ties are resolved by chance
df_s = df_s.drop_duplicates(subset='contig_CRISPR_spacer', keep='first')
# find unmatched spacers
spacers_with_match = set(df_s['contig_CRISPR_spacer'])
unmatched_spacers = all_spacer_names.difference(spacers_with_match)
# add unmatched spacers to the end of the spacer dataframe
for unmatched in unmatched_spacers:
    df_s = df_s.append(pd.Series({'contig_CRISPR_spacer': unmatched, 
                                  'accession': "no_match", 
                                  'common_name': None, 
                                  'e_value': None,
                                  'score': None,
                                  'percent_ident': None}), ignore_index=True)

# parse the contig, CRISPR, and spacer indices from the first column and then drop that column
df_s['contig_id'] = df_s['contig_CRISPR_spacer'].apply(lambda x: '_'.join(x.split('_')[:2]))
df_s['CRISPR'] = df_s['contig_CRISPR_spacer'].apply(lambda x: float(x.split('_')[2]))
df_s['spacer'] = df_s['contig_CRISPR_spacer'].apply(lambda x: float(x.split('_')[3]))
df_s = df_s.drop('contig_CRISPR_spacer', axis=1)

# merge with the contig dataframe, duplicating contig rows when they contained multiple spacers
df_out = df.merge(df_s, how='outer', on='contig_id', suffixes=('', '_spacer'))

# parse the second part of the contig name for a numerical sort (aesthetic)
df_out['contig_no'] = df_out['contig_id'].apply(lambda x: float(x.split('_')[1]))
df_out = df_out.sort_values(['contig_no', 'CRISPR', 'spacer'])
df_out = df_out.drop(['contig_no'], axis=1).set_index('contig_id')

# use pattern matching in the CCF Cas file to collect any Cas protein sequences for each contig
# place them in a temporary file
#with open(args.cas_file, 'r') as f:
#    with open(os.path.join(args.out_path, 'cas_list.txt'), 'w') as t:
#        for line in f.readlines():
#            if len(line.split('_')) > 2 and line[0] != '#':
#                contig = '_'.join(line.split('\t')[0].split('_')[:3]) 
#                cas = line.split('\t')[1]
#                t.write(f"{contig}\t{cas}\n")

# read this new cas file into a pandas dataframe
#df_cas = pd.read_csv(os.path.join(args.out_path, 'cas_list.txt'), delimiter='\t', header=None)
#df_cas.columns = ['contig_id', 'cas']

# add a new column which will allow us to see which types of Cas were present after reshaping
#df_cas['present'] = 1

# what's going on here? should be many duplicate contig_ids with different cas
#df_cas = df_cas.sort_values('contig_id').reset_index(drop=True).drop_duplicates(subset='contig_id')

#df_cas = df_cas.pivot(index='contig_id', columns='cas', values='present').reset_index()
#df_cas['contig_id'] = df_cas['contig_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
#df_cas = df_cas.fillna(0)
#df_cas = df_cas.groupby('contig_id').agg('sum')#.replace(0, None)
#df_out = df_out.merge(df_cas, how='left', on='contig_id')

df_out.to_csv(os.path.join(args.out_path, 'final_output.tab'), sep='\t')
