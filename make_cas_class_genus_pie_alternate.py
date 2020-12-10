import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(description='Create a 3-layer pie chart of detected cas genes, their classes, and their parent taxa')
parser.add_argument('--table', action='store', dest='table_path', help='full path to the final_output.tab')
parser.add_argument('--cas', action='store', dest='cas_path', help='full path to Cas_summary.tsv')
parser.add_argument('--out', action='store', dest='out_path')
args = parser.parse_args()

df = pd.read_csv(args.table_path, sep='\t', index_col=0)
df_cas = pd.read_csv(args.cas_path, sep='\t')
df = df.reset_index().drop_duplicates(subset='contig_id', keep='first')
df = df.merge(df_cas, on='contig_id', how='outer')

df['common_name'] = df['common_name'].fillna('Unknown unknown')
df['genus'] = df['common_name'].apply(lambda x: x.split(' ')[0].strip('[').strip(']'))
df['species'] = df['common_name'].apply(lambda x: x.split(' ')[1])
df = df[['contig_id', 'genus', 'species', 'cluster_type', 'cas_gene']]
df.loc[~df['cas_gene'].isna(), 'cas_count'] = 1
df['cas_count'] = df['cas_count'].fillna(0)
df = df.loc[df['cas_count']>0]
df['cas_gene'] = df['cas_gene'].apply(lambda x: x.split('_')[0])
df=df.reset_index()

df = df.groupby(['genus', 'contig_id', 'cas_gene', 'cluster_type'])['cas_count'].agg(['count']).reset_index()
df2 = df.groupby(['genus', 'contig_id'])['count'].sum().reset_index()
df3 = df.groupby(['genus'])['count'].sum().reset_index()

df = df.merge(df2, on=['genus', 'contig_id']).merge(df3, on='genus')
df.to_csv(os.path.join(args.out_path, 'cas_class_genus.tsv'))

fig, ax = plt.subplots()
cmap = plt.get_cmap('gist_rainbow')
cmap2 = plt.get_cmap('binary')

size = 0.4
wedges, _ = ax.pie(df.groupby('genus', sort=False).first()['count'].values, 
       radius=1-size, 
       wedgeprops=dict(width=size, edgecolor='w'), 
       colors=cmap(20*np.arange(len(df.groupby('genus', sort=False).first()['count'].values))),
       labels=list(df.groupby('genus', sort=False).first().index+': '+
                   df.groupby('genus', sort=False).first()['count'].astype(str)),
       rotatelabels=True,
       labeldistance=0.5,
       counterclock=False,
       startangle=-90)
ax.pie(list(df.groupby(['genus', 'contig_id'], sort=False).first()['count_y'].values), 
       radius=1, 
       wedgeprops=dict(width=size, edgecolor='w'), 
       colors=cmap2([100]*len(df.groupby(['genus', 'contig_id'], sort=False).first())),
       labels=list(df.groupby(['genus', 'contig_id']).first().reset_index()['contig_id']+': '+
                   df.groupby(['genus', 'contig_id']).first().reset_index()['count_y'].astype(str)),
                   rotatelabels=True,
                   labeldistance=0.65,
                   counterclock=False,
                   startangle=-90)
ax.pie(list(df['count_x'].values), 
       radius=1+size, 
       wedgeprops=dict(width=size, edgecolor='w'), 
       colors=cmap((df['cluster_type']=='General-Class1').astype(int)*100),
       labels=list(df['cas_gene']+': '+df['count_x'].astype(str)),
       rotatelabels=True,
       labeldistance=0.8,
       counterclock=False,
       startangle=-90)

ax.legend(wedges, df['cluster_type'].unique())
          
leg = ax.get_legend()
leg.legendHandles[0].set_color('lime')
leg.legendHandles[1].set_color('red')

figure= plt.gcf()
figure.set_size_inches(15, 15)
plt.savefig(os.path.join(args.out_path, 'cas_class_genus_alternate.png'))
