import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(description='Create a 2-layer pie chart of genus and species for blast hits on contigs')
parser.add_argument('--table', action='store', dest='table_path', help='full path to the final_output.tab')
parser.add_argument('--out', action='store', dest='out_path')
args = parser.parse_args()

df = pd.read_csv(args.table_path, sep='\t', index_col=0)
df = df.reset_index().drop_duplicates(subset='contig_id', keep='first')
df['common_name'] = df['common_name'].fillna('Unknown unknown')
df['genus'] = df['common_name'].apply(lambda x: x.split(' ')[0].strip('[').strip(']'))
df['species'] = df['common_name'].apply(lambda x: x.split(' ')[1])
df = df[['genus', 'species']]
df['species_count'] = 1

df_species = df.groupby(['genus', 'species']).count().reset_index()
df_out_1 = df_species[['genus', 'species', 'species_count']].sort_values('species_count', ascending=False).set_index(['genus', 'species'])
df_out_1.to_csv(os.path.join(args.out_path, 'species_contig_counts.tsv'), sep='\t')

fig, ax = plt.subplots()
cmap = plt.get_cmap('gist_rainbow')

size = 0.4
species_sum = df_species['species_count'].sum()
df_inner = df_species.groupby('genus').sum(numeric_only=True).sort_values(by='species_count', ascending=False).head(10)
df_inner.loc['Other', 'species_count'] = species_sum - sum(df_inner['species_count'])
df_inner = df_inner.rename(columns={'species_count': 'genus_count'})

df_outer = df_species.loc[df_species['genus'].isin(df_inner.index)]
df_outer = df_outer.append(pd.Series({'genus': 'Other','species': '','species_count': species_sum - sum(df_outer['species_count'])}),ignore_index=True)
df_outer = df_outer.sort_values(by='species_count', ascending=False)#.groupby('genus').head(5)
df_outer = df_outer.loc[df_outer['species_count']>2]
df_new = df_outer.merge(df_inner.reset_index(), on='genus', how='outer')
df_new['sort_key_outer'] = 1
df_new.loc[df_new['genus']=='Other', 'sort_key_outer'] = 0

other_species = df_new.groupby('genus')['genus_count'].first() - df_new.groupby('genus')['species_count'].sum()
df_new['sort_key'] = 1
for genus, val in other_species.to_dict().items():
    if val > 0:
        df_new = df_new.append(pd.Series({'genus': genus, 
                                          'species': 'other', 
                                          'species_count': val, 
                                          'sort_key': 0, 
                                          'sort_key_outer': 1,
                                          'genus_count': df_new.loc[df_new['genus']==genus, 'genus_count'].values[0]}), 
                               ignore_index=True)
df_new = df_new[['genus', 'genus_count', 'species', 'species_count', 'sort_key', 'sort_key_outer']].sort_values(['sort_key_outer', 'genus_count', 'sort_key', 'species_count'], ascending=False)
df_new[['genus', 'genus_count', 'species', 'species_count']].to_csv(os.path.join(args.out_path, 'taxonomy_breakdown.tsv'), sep='\t')

ax.pie(df_new['species_count'].values, 
    radius=1, 
    wedgeprops=dict(width=size, edgecolor='w'), 
    colors=cmap(7*np.arange(len(df_new['species_count'].values))),
    labels=list(df_new['species'] + ': ' + df_new['species_count'].astype(int).astype(str)),
    rotatelabels=True,
    labeldistance=1.05,
    explode= [0.1]*len(df_new['species_count']),
    counterclock=False,
    startangle=-90)
ax.pie(list(df_new.groupby('genus', sort=False).first()['genus_count'].values), 
    radius=1-size, 
    wedgeprops=dict(width=size, edgecolor='w'), 
    colors=cmap(27*np.arange(len(df_new.groupby('genus').first()['genus_count'].values))),
    labels=list(df_new.groupby('genus', sort=False).first().index + ': ' +
              ((df_new.groupby('genus', sort=False).first()['genus_count']/
                df_new.groupby('genus', sort=False).first()['genus_count'].sum())*100).round(1).astype(str) + '%'),
    rotatelabels=True,
    labeldistance=0.2,
    counterclock=False,
    startangle=-90)

figure= plt.gcf()
figure.set_size_inches(15, 15)
plt.savefig(os.path.join(args.out_path, 'genus_species_pie.png'))

