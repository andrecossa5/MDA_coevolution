"""
WES mutations for each clone.
"""

import os
import numpy as np
import pandas as pd
from Cellula.plotting._plotting_base import *


##


# Paths
path_main = '/Users/IEO5505/Desktop/MDA_coevolution_project'
path_data = os.path.join(path_main, 'data', 'WES')
path_results = os.path.join(path_main, 'results', 'WES')

# Read and format
L_all = []
L = []
for x in os.listdir(path_data):
    df_sample = pd.read_csv(os.path.join(path_data, x, f'filtered.{x}.small_mutations.cancervar.escat.maf.tsv'), sep='\t')
    df_sample = df_sample.rename(columns={'Tumor_Sample_Barcode':'clone'})
    L_all.append(df_sample)
    df_ = (
        df_sample
        .query('tumor_f>.1 and CancerVar == "Tier_II_potential" or ClinVar_VCF_CLNSIG == "Pathogenic"')
        [['Hugo_Symbol', 'Variant_Type', 'Variant_Classification', 'Variant_Type', 'cDNA_Change', 'tumor_f']]
    )
    df_.to_csv(os.path.join(path_results, f'{x}_vaf10_pathogenic.csv'))
    L.append(df_.assign(clone=x))

df_all = pd.concat(L_all)
df = pd.concat(L)


##


# VAF
fig, ax = plt.subplots(figsize=(5,5))
sns.kdeplot(df_all, x='tumor_f', hue='clone', fill=True, ax=ax)
ax.set(title='VAF distributions')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'VAF_distributions.png'), dpi=300)

# n drivers
fig, ax = plt.subplots(figsize=(5,5))
bar(df.groupby('clone').size().to_frame('n'), 'n', ax=ax, c='k', s=.75)
format_ax(ax, title='n putative "drivers"', xlabel='Clone', ylabel='n')
ax.spines[['right', 'top']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'n_drivers.png'), dpi=300)

# n common muts
clones = df['clone'].unique()
n = clones.size
M = np.zeros((n,n))

for i, x in enumerate(clones):
    df_x = df.query('clone == @x')
    muts_x = df_x['Hugo_Symbol'] + '_' + df_x['cDNA_Change']
    for j, y in enumerate(clones):
        df_y = df.query('clone == @y')
        muts_y = df_y['Hugo_Symbol'] + '_' + df_y['cDNA_Change']
        M[i,j] = len(set(muts_y.to_list()) & set(muts_x.to_list()))

fig, ax = plt.subplots(figsize=(5,5))
plot_heatmap(pd.DataFrame(M, index=clones, columns=clones), ax=ax, 
            title=f'Common drivers (total 13)', annot=True, y_names=False, label='n muts', annot_size=7)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'comon_drivers.png'), dpi=300)

