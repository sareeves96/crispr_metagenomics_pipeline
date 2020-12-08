from Bio import SeqIO, SearchIO
import argparse
import pandas as pd
import json
import os
import matplotlib.pyplot as plt
# Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table

parser = argparse.ArgumentParser(description='Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table')
parser.add_argument('--contigs', action='store', dest='contig_file_path', help='full path to original contigs fasta file')
parser.add_argument('--spacers', action='store', dest='spacer_file_path', help='full path to labelled_spacers.fa file')
parser.add_argument('--contigs_blast', action='store', dest='blast_result', help='full path to blast results tab from contig sequences')
parser.add_argument('--spacers_blast', action='store', dest='spacer_result', help='full path to blast results tab from spacer sequences')
parser.add_argument('--result', action='store', dest='json_file_path', help='full path to result.json in CCF output')
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

df_blast = df.copy(deep=True)

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

df_out.to_csv(os.path.join(args.out_path, 'final_output.tab'), sep='\t')

# Parse CRISPRCasFinder results into a Python dictionary
json_dict = []
with open(args.json_file_path, 'r') as file:
    json_dict = json.load(file)
    file.close()

# Only useful data is contained in sequences
json_seq_list = json_dict['Sequences']

#Write a row for each cas gene detected, associated with a type and contig_id
with open(os.path.join(args.out_path, 'Cas_summary.tsv'), 'w') as f:
    f.write(f"contig_id\tcluster_type\tcas_gene\n")
    for seq in json_seq_list:
        id_ = seq["Id"]
        for cas_cluster in seq['Cas']:
            cas_type = cas_cluster['Type']
            for gene in cas_cluster['Genes']:
                gene_subtype = gene['Sub_type']
                f.write(f"{id_}\t{cas_type}\t{gene_subtype}\n")

# read this new cas file into a pandas dataframe
df_cas = pd.read_csv(os.path.join(args.out_path, 'Cas_summary.tsv'), delimiter='\t')
df_cas = df_cas.drop('cluster_type', axis=1)
df_cas['cas_gene'] = df_cas['cas_gene'].apply(lambda x: x.split('_')[0])
#df_cas.columns = ['contig_id', 'cas_gene']

# add a new column which will allow us to see which types of Cas were present after reshaping
df_cas['present'] = 1

df_cas = df_cas.sort_values('contig_id').reset_index(drop=True)
df_cas = df_cas.groupby(['contig_id', 'cas_gene']).sum().reset_index()
df_cas = df_cas.pivot(index='contig_id', columns='cas_gene', values='present').reset_index() # takes each cas_gene as a column and the contig id's as rows/indices and fills the cell with 1 if the cas_gene is present in the contig
df_cas['contig_id'] = df_cas['contig_id'].apply(lambda x: '_'.join(x.split('_')[:2])) # label formatting
df_cas = df_cas.fillna(0)  # fills empty cells with 0
df_cas = df_cas.groupby('contig_id').agg('sum')
#print(df_cas)

## add the number of occurences of different cas genes, with each column representing a cas gene, 1 x n pandas dataframe,
##     n being the number of unique cas genes identified.
df_sum_cas = df_cas.agg(['sum'])
#print(df_sum_cas)

## assign the value of the total cas genes identified, all inclusive, to total_cas_genes.
total_cas_genes = df_sum_cas.agg(['sum'], axis=1)
total_cas_genes = total_cas_genes.to_numpy()[0][0]   # total number of cas genes identifed in the metagenome (sum of the values in row/index 'sum' in df_sum_cas.
#print("total cas genes in metagenome: ", total_cas_genes)

## Because these are counts
for col in df_sum_cas.columns:
    df_sum_cas[col] = df_sum_cas[col].astype(int)

## calculate the occurence of each cas gene as a percentage of total cas genes identified in the metagenome
df_casgene_occurence = df_sum_cas.apply(lambda x: (x/total_cas_genes)*100).round(2).T  # rounds the percentages of cas genes to 2 d.p
df_casgene_occurence.columns = ['Occurence']
print(df_casgene_occurence)

plot_title = "Distribution of cas genes in the metagenome as a percentage"
plot_pie = df_casgene_occurence.plot.pie(y='Occurence', figsize=(8,8), ylabel="", legend=False, autopct='%.2f', fontsize=10, title=plot_title)  # Plot
pie_chart = plot_pie.get_figure()
plt.tight_layout()
pie_chart.savefig(os.path.join(args.out_path, 'cas_pie_chart.pdf'))  # exports plot

## --------------------------------------------------

## --------------------------------------------------
## Export the frequency of occurence of all cas genes as a matplotlib table to a HTML file and as a tab seperated file

df_cas_transposed = df_sum_cas.transpose()
df_cas_transposed.columns = ['Counts']
df_cas_transposed.to_csv(os.path.join(args.out_path, 'casgene_occurence.tab'), sep='\t')

## **************************************************--------------------------------------------------
## **************************************************
## Distribution of cas genes with respect to the genera in the metagenome
## --------------------------------------------------

## Clean up
#df_contigs_cas = df_cas.agg(['sum'], axis=1).reset_index()

df = df_out
df['common_name'] = df['common_name'].fillna('Unknown unknown')
df['genus'] = df['common_name'].apply(lambda x: x.split(' ')[0].strip('[').strip(']'))
df['species'] = df['common_name'].apply(lambda x: x.split(' ')[1])
df_genera = df[['genus', 'species']]

df_cas_genera = df_cas.join(df_genera).reset_index().drop_duplicates(subset='contig_id').set_index('contig_id')
print(len(df_cas_genera))
#print(df_genera_cas)
#print(len(df_genera_cas))

#df_cas_genera = pd.DataFrame.from_dict(dict_cas_genera, orient='index', columns=['total_freq']) # genera are the indices
df_cas_genera = df_cas_genera.reset_index().set_index(['contig_id', 'genus', 'species'])
df_cas_genera['total'] = df_cas_genera.values.sum(axis=1)

total_cas_genera = df_cas_genera['total'].sum()
## Occurence of cas genes in various genera as a percentage of total cas gene occurences in the metagenome
df_genus_cas_occurence = df_cas_genera['total'].apply(lambda x: (x/total_cas_genera)*100).round(2).reset_index()

df_p= df_genus_cas_occurence.groupby('genus').sum()
df_pie = df_p.loc[df_p['total']>2]
df_pie.loc['Other', 'total'] = df_p.loc[df_p['total']<=2].sum().values[0]
df_pie = df_pie.rename(columns={'total':'Occurence'})
print(df_pie)

plot_title_cas_genera = "Distribution of cas genes across genera as a percentage"
plot_pie = df_pie.plot.pie(y='Occurence', figsize=(8,8), ylabel="", legend=False, autopct='%.2f', fontsize=10, title=plot_title_cas_genera)  # Plot
pie_chart_cas_genus = plot_pie.get_figure()
plt.tight_layout()
pie_chart_cas_genus.savefig(os.path.join(args.out_path, 'genera_with_cas_pie_chart.pdf'))  # exports plot
df_pie.to_csv(os.path.join(args.out_path, 'cas_genera_distribution.tab'), sep='\t', index_label='Genus')

## **************************************************--------------------------------------------------
## **************************************************
## Ratio of the number of species with at least 1 cas gene in a genus and the total identified species in that genus in the metagenome
## --------------------------------------------------
df_outer = df_genera.join(df_cas).reset_index().drop_duplicates(subset='contig_id').set_index(['contig_id', 'genus', 'species']).fillna(0)
df_outer['total'] = df_outer.values.sum(axis=1)
df_outer['total'] = df_outer['total'].apply(lambda x: 1 if x>0 else 0)
df_outer = df_outer['total'].reset_index().groupby('genus')['total'].agg(['count', 'sum'])
df_ratios = df_outer.loc[(df_outer['count']>1)&(df_outer['sum']>0)]
df_ratios['percent'] = df_ratios['sum']/df_ratios['count']*100
df_ratios = df_ratios[['percent']].reset_index()
## --------------------------------------------------
## Plot the data as a pie chart and export as a PDF
plt.clf()
plot_title_cas_genera = "Fractions of genera with cas genes"
plot_bar = df_ratios.plot.bar(y='percent', x='genus', figsize=(8,5), ylabel="Ratios in percentage", rot=90, legend=False, fontsize=10, title=plot_title_cas_genera)  # Plot
bar_chart_ratios = plot_bar.get_figure()
plt.tight_layout()
bar_chart_ratios.savefig(os.path.join(args.out_path, 'ratios_bar_chart.pdf'), orientation='landscape')  # exports plot

## --------------------------------------------------
## Export the ratios as a table to a HTML file and as a tab seperated file

df_ratios.to_csv(os.path.join(args.out_path, 'cas_occurence_genera.tab'), sep='\t', index_label='Genus')
