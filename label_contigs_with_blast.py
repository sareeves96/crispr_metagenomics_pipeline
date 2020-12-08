from Bio import SeqIO, SearchIO
import argparse
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from xhtml2pdf import pisa
from PyPDF4 import PdfFileMerger


# Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table

parser = argparse.ArgumentParser(description='Connects contig identities, taxonomies, spacers, and spacer taxonomies into one table')
parser.add_argument('--contigs', action='store', dest='contig_file_path', help='full path to original contigs fasta file')
parser.add_argument('--spacers', action='store', dest='spacer_file_path', help='full path to labelled_spacers.fa file')
parser.add_argument('--contigs_blast', action='store', dest='blast_result', help='full path to blast results tab from contig sequences')
parser.add_argument('--spacers_blast', action='store', dest='spacer_result', help='full path to blast results tab from spacer sequences')
parser.add_argument('--cas_file', action='store', dest='cas_file', help='full path to Cas.REPORT.tsv file in the CCF output directory')
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




## **************************************************--------------------------------------------------
## **************************************************
## Distribution of cas genes in the metagenome
## --------------------------------------------------

## use pattern matching in the CCF Cas file (Cas.REPORT.tsv) to collect any identified cas gene sequences for each contig 
## place them in a temporary file

with open(args.cas_file, 'r') as f:
    with open(os.path.join(args.out_path, 'temp/temp_cas_genes.txt'), 'w') as t:
        for line in f.readlines():
            if len(line.split('_')) > 2 and line[0] != '#':
                contig = '_'.join(line.split('\t')[0].split('_')[:3]) 
                cas = line.split('\t')[1].split("_") # the Cas gene labels are in the following format: cas1_0_IE, cas2_0_I-II-III, cas2_0_I-II-III-IV, cas2_0_IE, etc..
                cas_gene = cas[0]
                t.write(f"{contig}\t{cas_gene}\n")

# read this new cas file into a pandas dataframe
df_cas = pd.read_csv(os.path.join(args.out_path, 'temp/temp_cas_genes.txt'), delimiter='\t', header=None)
df_cas.columns = ['contig_id', 'cas_gene']

# add a new column which will allow us to see which types of Cas were present after reshaping
df_cas['present'] = 1

df_cas = df_cas.sort_values('contig_id').reset_index(drop=True) 
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
print("total cas genes in metagenome: ", total_cas_genes)

## Because these are counts
for col in df_sum_cas.columns:
    df_sum_cas[col] = df_sum_cas[col].astype(int)

## calculate the occurence of each cas gene as a percentage of total cas genes identified in the metagenome
df_casgene_occurence = df_sum_cas.apply(lambda x: (x/total_cas_genes)*100).round(2)  # rounds the percentages of cas genes to 2 d.p

#print(df_casgene_occurence)
#print(df_casgene_occurence.index)
#print(df_casgene_occurence.columns)

## --------------------------------------------------
## highlight the non-trivial data 

def condense(dict_name, in_df_name, col1, col2, threshold, ind, col_row):
    '''
    Creates a python dictionary (dict_name) with two keys (col1 and col2). Assign the labels to col1 (will be converted to indeces later)
        and assign the occurence to col 2. Set 'ind' to True if iterating over rows/indices in the for loop, False otherwise. col_row is the
        non-iterating column or row.
    
    condense: Str (pandas DataFrame) Str Str (anyof Int Float) Bool Str -> dict
    '''
    others = 0
    dict_name = {}
    dict_name[col1] = []
    dict_name[col2] = []
    if ind:
        for row in in_df_name.index:
            occurence = in_df_name.loc[row, col_row]
            if occurence < threshold:
                others += occurence
            else:
                dict_name[col1].append(row)
                dict_name[col2].append(occurence)
    else:
        for col in in_df_name.columns:
            occurence = in_df_name.loc[col_row, col]
            if occurence < threshold:
                others += occurence
            else:
                dict_name[col1].append(col)
                dict_name[col2].append(occurence)

    dict_name[col1].append('others')
    dict_name[col2].append(others)
    return dict_name

## --------------------------------------------------

dict_pie_cas = condense("dict_pie_cas", df_casgene_occurence, "cas_gene", "Occurence", 5, False, "sum")
df_pie_cas = pd.DataFrame.from_dict(dict_pie_cas) 

## --------------------------------------------------
## Plot the data as a pie chart and exports it as a PDF file

df_pie_cas = pd.DataFrame.from_dict(dict_pie_cas)
df_pie_cas.pop('cas_gene')
df_pie_cas.index = dict_pie_cas['cas_gene']
print(df_pie_cas)

plot_title = "Distribution of cas genes in the metagenome as a percentage"
plot_pie = df_pie_cas.plot.pie(y='Occurence', figsize=(8,8), ylabel="", legend=False, autopct='%.2f', fontsize=10, title=plot_title)  # Plot
pie_chart = plot_pie.get_figure()
plt.tight_layout()
pie_chart.savefig(os.path.join(args.out_path, 'cas_pie_chart.pdf'))  # exports plot

## --------------------------------------------------

## --------------------------------------------------
## Export the frequency of occurence of all cas genes as a matplotlib table to a HTML file and as a tab seperated file

df_cas_transposed = df_sum_cas.transpose()
df_cas_transposed.columns = ['Counts']
# print(df_cas_transposed)
df_cas_transposed.to_html(os.path.join(args.out_path, 'temp/temp_cas_table.html'))
df_cas_transposed.to_csv(os.path.join(args.out_path, 'casgene_occurence.tab'), sep='\t', index_label='cas_gene')

## --------------------------------------------------

## --------------------------------------------------
## Merge all data into 1 PDF file

def convert_to_pdf(file_name):
    '''
    Converts a HTML file to PDF. Returns None. Do not include an extension for file_name; extensions are added within the function.

    convert_to_pdf: Str -> None

    Requires: Input file is a HTML file
    '''
    html_file = open(os.path.join(args.out_path, (file_name + ".html")), 'r+b')
    result_file = open(os.path.join(args.out_path, (file_name + '.pdf')), 'w+b')
    pisa_status = pisa.CreatePDF(html_file, dest=result_file)
    result_file.close()
    html_file.close()

def merge_pdfs(pdfs, out):
    '''
    Merges all PDFs in pdfs and exports as one PDF file with out as the label. Specify extension (.pdf). Returns None.

    merge_pdfs: (listof Str) Str -> None

    Requires: all strings in pdfs are referring to a PDF file.
    '''
    all_pdfs = []
    for pdf in pdfs:
        read = open(os.path.join(args.out_path, pdf), 'r+b')
        all_pdfs.append(read)

    merger = PdfFileMerger()
    for pdf in all_pdfs:
        merger.append(pdf)

    merger.write(os.path.join(args.out_path, out))
    merger.close()

## --------------------------------------------------

## convert the exported HTML table to a PDF file.
convert_to_pdf('temp/temp_cas_table')

## merge the pie chart and the table into one PDF file
cas_pdfs = ['cas_pie_chart.pdf', 'temp/temp_cas_table.pdf']
merge_pdfs(cas_pdfs, "casgene_distribution.pdf")

## --------------------------------------------------
## **************************************************


## **************************************************--------------------------------------------------
## **************************************************
## Distribution of genera in the metagenome
## --------------------------------------------------

## clean it up
df_blast.pop('accession')
df_blast.pop('e_value')
df_blast.index = df_blast['contig_id']

## isolate the labelled genus and species names for each contig and insert them as columns back into df_blast
genera = []
species = []
for ind in df_blast.index:
    if df_blast['common_name'][ind] == None:   # one of the contigs was unidentified in my data
        g = "unidentified"
        s = "species"
    else:
        sciname = df_blast['common_name'][ind].split(" ")
        g = sciname[0].strip('[]')
        s = sciname[1]
    genera.append(g)
    species.append(s)

df_blast.insert(loc=1, column='Genus', value=genera)
df_blast.insert(loc=2, column='species', value=species)

species2 = species[:]

## Count occurence of each genus.
## Accounted for multiple contigs originating from the same species, assuming blast labelling was accurate. 
dict_genera = {}
for ind in df_blast.index:
    genus = df_blast.loc[ind, 'Genus']
    sp = df_blast.loc[ind, 'species']
    if genus in dict_genera:
        if sp in species:
            dict_genera[genus].append(1)
            species.remove(sp)
    else:
        if sp in species:
            dict_genera[genus] = []
            dict_genera[genus].append(1)
            species.remove(sp)

for key in dict_genera.keys():
    dict_genera[key] = sum(dict_genera[key])

## Create a pandas df using the dictionary
df_genera = pd.DataFrame.from_dict(dict_genera, orient='index', columns=['total_freq']) # orient='index' changes the dictionary columns to rows in the df

## total number of different genera identified
total_genera = df_genera.agg(['sum'])['total_freq']['sum']

## Occurence of each genus in the metagenome
df_genus_occurence = df_genera.apply(lambda x: (x/total_genera)*100).round(2)

## Highlight non-trivial data. Threshold = 1
dict_pie_genera = condense("dict_pie_genera", df_genus_occurence, "Genus", "Occurence", 1, True, "total_freq")

## --------------------------------------------------
## Plot the data as a pie chart and exports it as a PDF file

df_pie_genera = pd.DataFrame.from_dict(dict_pie_genera)

plot_title_genera = "Distribution of genera in the metagenome as a percentage"
plot_pie = df_pie_genera.plot.pie(y='Occurence', x='Genus', figsize=(8,8), ylabel="", legend=False, autopct='%.2f', fontsize=10, title=plot_title_genera)  # Plot
pie_chart_genera = plot_pie.get_figure()
plt.tight_layout()
pie_chart_genera.savefig(os.path.join(args.out_path, 'genera_pie_chart.pdf'))  # exports plot

## --------------------------------------------------


## --------------------------------------------------
## Export the frequency of occurence of all genera as a matplotlib table to a HTML file and as a tab seperated file

df_genera.to_html(os.path.join(args.out_path, 'temp/temp_all_genera_table.html'))
df_genera.to_csv(os.path.join(args.out_path, 'genera_occurence.tab'), sep='\t', index_label='cas_gene')

## Covert the HTML file to PDF
convert_to_pdf('temp/temp_all_genera_table')
## --------------------------------------------------
## **************************************************




## **************************************************--------------------------------------------------
## **************************************************
## Distribution of cas genes with respect to the genera in the metagenome
## --------------------------------------------------

## Clean up
df_contigs_cas = df_cas.agg(['sum'], axis=1)

## isolate contig ids, genus and species names of contigs with identified cas genes 
dict_genera_cas = {}
dict_genera_cas['contig_id'] = []
dict_genera_cas['Genus'] = []
dict_genera_cas['species'] = []
for ind in df_contigs_cas.index:
    for ind_blast in df_blast.index:
        if ind == ind_blast:
            dict_genera_cas['contig_id'].append(ind)
            dict_genera_cas['Genus'].append(df_blast.loc[ind_blast, 'Genus'])
            dict_genera_cas['species'].append(df_blast.loc[ind_blast, 'species'])

df_genera_cas = pd.DataFrame.from_dict(dict_genera_cas)

## Count the occurence of cas genes within each genus
dict_cas_genera = {}
for ind in df_genera_cas.index:
    genus = df_genera_cas.loc[ind, 'Genus']
    sp = df_genera_cas.loc[ind, 'species']
    if genus in dict_cas_genera.keys():
        if sp in species2:
            dict_cas_genera[genus].append(1)
            species2.remove(sp)
    else:
        if sp in species2:
            dict_cas_genera[genus] = []
            dict_cas_genera[genus].append(1)
            species2.remove(sp)

for key in dict_cas_genera.keys():
    dict_cas_genera[key] = sum(dict_cas_genera[key])

df_cas_genera = pd.DataFrame.from_dict(dict_cas_genera, orient='index', columns=['total_freq']) # genera are the indices

total_cas_genera = df_cas_genera.agg(['sum'])['total_freq']['sum']

## Occurence of cas genes in various genera as a percentage of total cas gene occurences in the metagenome
df_genus_cas_occurence = df_cas_genera.apply(lambda x: (x/total_cas_genera)*100).round(2)

#print(df_genus_cas_occurence)


## Highlight non-trivial data. Threshold = 2%
dict_pie_cas_genus = condense("dict_pie_cas_genus", df_genus_cas_occurence, "Genus", "Occurence", 2, True, "total_freq")

## --------------------------------------------------
## Plot the data as a pie chart and export as a PDF

df_pie_cas_genus = pd.DataFrame.from_dict(dict_pie_cas_genus)
df_pie_cas_genus.index = df_pie_cas_genus['Genus']


plot_title_cas_genera = "Distribution of cas genes across genera as a percentage"
plot_pie = df_pie_cas_genus.plot.pie(y='Occurence', x='Genus', figsize=(8,8), ylabel="", legend=False, autopct='%.2f', fontsize=10, title=plot_title_cas_genera)  # Plot
pie_chart_cas_genus = plot_pie.get_figure()
plt.tight_layout()
pie_chart_cas_genus.savefig(os.path.join(args.out_path, 'genera_with_cas_pie_chart.pdf'))  # exports plot

## --------------------------------------------------


## --------------------------------------------------
## Export the frequency of occurence of all cas genes with respect the genera as a matplotlib table to a HTML file and as a tab seperated file

df_cas_genera.to_html(os.path.join(args.out_path, 'temp/temp_cas_genera_table.html'))
df_cas_genera.to_csv(os.path.join(args.out_path, 'genera_cas_occurence.tab'), sep='\t', index_label='Genus')

## Covert the HTML file to PDF
convert_to_pdf('temp/temp_cas_genera_table')
## --------------------------------------------------
## **************************************************




## **************************************************--------------------------------------------------
## **************************************************
## Ratio of the number of species with at least 1 cas gene in a genus and the total identified species in that genus in the metagenome
## --------------------------------------------------

## Calculate the ratios in percentage
dict_calc_ratios = {}
dict_calc_ratios['Genus'] = []
dict_calc_ratios['Ratio'] = []
for cas_genus in df_cas_genera.index:
    cas_genus_count = df_cas_genera.loc[cas_genus, 'total_freq']
    for genus in df_genera.index:
        if cas_genus == genus:
            genus_count = df_genera.loc[genus, 'total_freq']
            ratio_genus = (cas_genus_count/genus_count) * 100
            dict_calc_ratios['Genus'].append(cas_genus)
            dict_calc_ratios['Ratio'].append(ratio_genus)


df_genus_ratios = pd.DataFrame.from_dict(dict_calc_ratios).round(2)
df_genus_ratios.index = dict_calc_ratios['Genus']
df_genus_ratios.pop('Genus')

df_ratios = df_genus_ratios.transpose()
df_ratios.index = ['Genera']

## --------------------------------------------------
## Plot the data as a pie chart and export as a PDF

plot_title_cas_genera = "Fractions of genera with cas genes"
plot_bar = df_ratios.plot.bar(figsize=(8,5), ylabel="Ratios in percentage", rot=0, legend=True, fontsize=10, title=plot_title_cas_genera, xticks=[])  # Plot
bar_chart_ratios = plot_bar.get_figure()
plt.tight_layout()
bar_chart_ratios.savefig(os.path.join(args.out_path, 'ratios_bar_chart.pdf'), orientation='landscape')  # exports plot

## --------------------------------------------------
## Export the ratios as a table to a HTML file and as a tab seperated file

df_genus_ratios.to_html(os.path.join(args.out_path, 'temp/temp_ratios_table.html'))
df_genus_ratios.to_csv(os.path.join(args.out_path, 'genera_cas_fraction.tab'), sep='\t', index_label='Genus')

## Convert HTML to PDF
convert_to_pdf('temp/temp_ratios_table')

## --------------------------------------------------
## --------------------------------------------------


## --------------------------------------------------
## Merge all the information concerning genera into one PDF file

genera_pdfs = ['genera_pie_chart.pdf', 'temp/temp_all_genera_table.pdf',  'genera_with_cas_pie_chart.pdf', 'temp/temp_cas_genera_table.pdf', 'ratios_bar_chart.pdf', 'temp/temp_ratios_table.pdf']
merge_pdfs(genera_pdfs, 'genera_stats.pdf')
        
## --------------------------------------------------
## **************************************************

## **************************************************--------------------------------------------------
