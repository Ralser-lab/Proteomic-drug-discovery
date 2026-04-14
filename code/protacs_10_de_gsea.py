#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_10_de_gsea.py
Description:
    Perform GSEA on limma differential expression results from the HBD drug
    proteome screen at 1 µM and 10 µM, summarize enriched GO terms, and generate 
    plots including a probability mass function plot (PMF), NES-colored strip plot, 
    and dose–response line plots across concentrations of a target GO term. 

Author: Shaon Basu
Date: 2025-09-19

Inputs
------
- data/Cluster{1..15}_10uM_Limma_250305a.csv
- data/Cluster{1..15}_1p0uM_Limma_250305a.csv
- data/Cluster_LFCxPval_0p1uM_250305a.csv
- data/Cluster_LFCxPval_1uM_250305a.csv
- data/Cluster_LFCxPval_10uM_250305a.csv
- data/c5.go.cc.v2023.2.Hs.symbols.gmt

Outputs
-------
- data/gsea_FDR_NES_redux_HBDlib_250205a.csv
- figures/protacs_10_de_gsea_pmf.pdf
- figures/protacs_10_de_gsea_chemicalseries_1and10uM.pdf
- figures/protacs_10_<gsea_term>_lineplot.pdf

Requirements
------------
Python >= 3.8
Dependencies: pandas, numpy, matplotlib, seaborn, gseapy

"""
# %% Import modules
import pandas as pd
import matplotlib.pyplot as plt
from device_supportfunctions import GBDTUtils
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import seaborn as sns
import numpy as np
import gseapy as gp
import os

# Set relative paths
dir_main = os.path.dirname(__file__)
data_dir = os.path.join(dir_main, '..', 'data')
fig_out = os.path.join(dir_main, '..', 'figures')
GBDTUtils.configure_font()

def run_gsea(x, set = 'c5.go.cc.v2023.2.Hs.symbols.gmt'):
    """
    Run preranked GSEA on a ranked gene list.

    Parameters
    ----------
    x : pandas.Series
        Ranked gene list (index = genes, values = ranking scores).

    Returns
    -------
    pandas.DataFrame
        GSEA results table.
    """  
    gsea_results = gp.prerank(rnk=x,
                            gene_sets=os.path.join(data_dir, set),
                            outdir=data_dir,
                            min_size=15,  # Minimum size of gene set to consider
                            max_size=500, # Maximum size of gene set to consider
                            processes=4)
    gsea_df = gsea_results.res2d
    return gsea_df

def gsealimma(x):
    """
    Rank genes from limma results and run GSEA.

    Parameters
    ----------
    x : pandas.DataFrame
        Limma DE results with columns 'logFC' and 'adj.P.Val'.

    Returns
    -------
    pandas.DataFrame
        GSEA results table.
    """

    x['sign'] = np.sign(x['logFC']) * -np.log10(x['adj.P.Val'])
    ranked = x['sign'].sort_values(ascending=False)
    #ranked = x['logFC'].sort_values(ascending=False)
    gsea_de = run_gsea(ranked)
    return gsea_de


# Load differential expression profiles by chemical series in the HBD drug proteome screen
# Contrasts = Chemical series vs DMSO [limma]
# Extract 10 micromolar values

de_1 = pd.read_csv(os.path.join(data_dir, 'Cluster1_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0, header = 0)
de_2 = pd.read_csv(os.path.join(data_dir, 'Cluster2_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_3 = pd.read_csv(os.path.join(data_dir, 'Cluster3_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_4 = pd.read_csv(os.path.join(data_dir, 'Cluster4_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_5 = pd.read_csv(os.path.join(data_dir, 'Cluster5_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_6 = pd.read_csv(os.path.join(data_dir, 'Cluster6_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_7 = pd.read_csv(os.path.join(data_dir, 'Cluster7_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0, header = 0)
de_8 = pd.read_csv(os.path.join(data_dir, 'Cluster8_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_9 = pd.read_csv(os.path.join(data_dir, 'Cluster9_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_10 = pd.read_csv(os.path.join(data_dir, 'Cluster10_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_11 = pd.read_csv(os.path.join(data_dir, 'Cluster11_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_12 = pd.read_csv(os.path.join(data_dir, 'Cluster12_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_13 = pd.read_csv(os.path.join(data_dir, 'Cluster13_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_14 = pd.read_csv(os.path.join(data_dir, 'Cluster14_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_15 = pd.read_csv(os.path.join(data_dir, 'Cluster15_10uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)

# Extract 1 micromolar values 
de_1_1 = pd.read_csv(os.path.join(data_dir, 'Cluster1_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0, header = 0)
de_1_2 = pd.read_csv(os.path.join(data_dir, 'Cluster2_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_3 = pd.read_csv(os.path.join(data_dir, 'Cluster3_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_4 = pd.read_csv(os.path.join(data_dir, 'Cluster4_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_5 = pd.read_csv(os.path.join(data_dir, 'Cluster5_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_6 = pd.read_csv(os.path.join(data_dir, 'Cluster6_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_7 = pd.read_csv(os.path.join(data_dir, 'Cluster7_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0, header = 0)
de_1_8 = pd.read_csv(os.path.join(data_dir, 'Cluster8_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_9 = pd.read_csv(os.path.join(data_dir, 'Cluster9_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_10 = pd.read_csv(os.path.join(data_dir, 'Cluster10_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_11 = pd.read_csv(os.path.join(data_dir, 'Cluster11_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_12 = pd.read_csv(os.path.join(data_dir, 'Cluster12_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_13 = pd.read_csv(os.path.join(data_dir, 'Cluster13_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_14 = pd.read_csv(os.path.join(data_dir, 'Cluster14_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)
de_1_15 = pd.read_csv(os.path.join(data_dir, 'Cluster15_1p0uM_Limma_250305a.csv'),
                    delimiter = ';', decimal = ',', index_col=0)

# %% Run GSEA on limma DE @ 1 and 10 micromolar 

# GSEA on limma differential expression sets
df1 = gsealimma(de_1)
df2 = gsealimma(de_2)
df3 = gsealimma(de_3)
df4 = gsealimma(de_4)
df5 = gsealimma(de_5)
df6 = gsealimma(de_6)
df7 = gsealimma(de_7)
df8 = gsealimma(de_8)
df9 = gsealimma(de_9)
df10 = gsealimma(de_10)
df11 = gsealimma(de_11)
df12 = gsealimma(de_12)
df13 = gsealimma(de_13)
df14 = gsealimma(de_14)
df15 = gsealimma(de_15)

# Limma DE 1 micromolar 
df16 = gsealimma(de_1_1)
df17 = gsealimma(de_1_2)
df18 = gsealimma(de_1_3)
df19 = gsealimma(de_1_4)
df20 = gsealimma(de_1_5)
df21 = gsealimma(de_1_6)
df22 = gsealimma(de_1_7)
df23 = gsealimma(de_1_8)
df24 = gsealimma(de_1_9)
df25 = gsealimma(de_1_10)
df26 = gsealimma(de_1_11)
df27 = gsealimma(de_1_12)
df28 = gsealimma(de_1_13)
df29 = gsealimma(de_1_14)
df30 = gsealimma(de_1_15)

# %% Plot GSEA by NES and FDR cutoffs using 10 micromolar dataset
gsea_dict = {
    "df16": df16,
    "df17": df17,
    "df18": df18,
    "df19": df19,
    "df20": df20,
    "df21": df21,
    "df22": df22,
    "df23": df23,
    "df24": df24,
    "df25": df25,
    "df26": df26,
    "df27": df27,
    "df28": df28,
    "df29": df29,
    "df30": df30,
    "df1": df1,
    "df2": df2,
    "df3": df3,
    "df4": df4,
    "df5": df5,
    "df6": df6,
    "df7": df7,
    "df8": df8,
    "df9": df9,
    "df10": df10,
    "df11": df11,
    "df12": df12,
    "df13": df13,
    "df14": df14,
    "df15": df15}

# Create dataframe(dogma): Df containing all the gsea pathways
# Set up global filtering variables to create the df
gsea_df = pd.DataFrame()
cat_order = list(gsea_dict.keys())
filter = 0.01 #FDR cutoff
filter2 = 0 #NES cutoff
top = 5

# Create top GSEA table
for idx, df in gsea_dict.items():
    df_subset = df.loc[df['FDR q-val'] < filter]
    df_subset = df_subset.loc[abs(df_subset['NES']) > filter2]
    df_subset.sort_values(by='FDR q-val', ascending=True)
    df_subset = df_subset.head(top)
    df_subset.index = df_subset['Term']
    df_subset['condition'] = idx
    gsea_df = pd.concat([gsea_df,df_subset], axis =0)
gsea_term_df = gsea_df.copy()

# Split and take substring
gsea_term_df['Term'] = [idx.split('_', 1)[1] if '_' in idx else idx for idx in gsea_term_df['Term']]

# Term filter
gsea_term_df['Term'] = ['_'.join(idx.split('_')[-9:]) if idx.count('_') >= 6 else idx for idx in gsea_term_df['Term']]

# Lower case 
gsea_term_df['Term'] = gsea_term_df['Term'].str.lower().str.replace('_', ' ')

# Convert NES to numeric
gsea_term_df['NES'] = pd.to_numeric(gsea_term_df['NES'], errors = 'coerce')

# Sort df by NES and save
gsea_term_df.sort_values(by = 'NES', inplace = True, ascending = False)
gsea_term_df.to_csv(os.path.join(data_dir, 'gsea_FDR_NES_redux_HBDlib_250205a.csv'))

# Import gsea without doing above, and format for PMF
gsealist = pd.read_csv(os.path.join(data_dir,'gsea_FDR_NES_redux_HBDlib_250205a.csv'), index_col = 0)

# Create df with counts, proba and nes
pmf = pd.DataFrame({'count': gsealist['Term.1'].value_counts(),
                    'nes': gsealist.groupby(['Term.1'])['NES'].mean()})

pmf['proba'] = pmf['count']/pmf['count'].sum()

pmf.sort_values(['proba'], inplace = True, ascending = False)

# Plot pmf for gsea terms
fig, ax = plt.subplots(figsize =[6,5])
plt.title('Degrader library \nproteomic signatures', fontsize=16)
sns.barplot(x = pmf.index, y = pmf['proba'], hue = pmf['nes'], palette = 'RdBu_r')
ax.tick_params(axis='x', labelsize=15, rotation=90)  # Adjust font size and rotation
ax.tick_params(axis='y', labelsize=15)  # Adjust font size
ax.set_xlabel('')  # This removes the x-axis label
ax.set_ylabel('')
ax.legend(loc='upper right', bbox_to_anchor=(1,1), fontsize=8, title='NES', title_fontsize=10, handlelength=0.5, labelspacing=0.5, borderpad=0.5)
plt.savefig(os.path.join(fig_out, 'protacs_10_de_gsea_pmf.pdf'))
# plt.show()
plt.close()

# %% Plot the results as dot plot
norm = TwoSlopeNorm(vmin=gsealist['NES'].min(), vcenter=0, vmax=gsealist['NES'].max())
gradient = 'RdBu_r'
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'

# Cluster name labels (series 1-15, ordered by PC1 distance)
cluster_labels = {
    1: 'AR-VHL-Other',    2: 'AR-T6N-Indole',    3: 'AR-Other-Other',
    4: 'Txn-VHL/6N-Other', 5: 'AR-VHL-Indole',   6: 'AR-T5N-Other',
    7: 'Txn-Other-Other',  8: 'AR-T6N/VHL-Pip',  9: 'AR-T6N-Other',
    10: 'AR-T5N-Pip',     11: 'AR-DHU-Pip',       12: 'AR-Other-Pip',
    13: 'Non PROTAC',     14: 'AR-L5N-Other',     15: 'AR-L5N-Pip'
}
# df16-df30 = clusters 1-15 at 1µM; df1-df15 = clusters 1-15 at 10µM
label_map = {f'df{i+16}': f'{cluster_labels[i+1]}\n(1µM)' for i in range(15)}
label_map.update({f'df{i+1}': f'{cluster_labels[i+1]}\n(10µM)' for i in range(15)})

# Create a stripplot (dot plot)
gseaplot = plt.figure(figsize=(30, 10))
ax = plt.gca()
sns.stripplot(x="condition", y="Term.1", data=gsealist, hue="NES", order=cat_order,
            palette=gradient, jitter=False, size=30, dodge=False, ax=ax)
# Recolor dots using TwoSlopeNorm
cmap_obj = plt.cm.get_cmap(gradient)
for coll_idx, cond in enumerate(cat_order):
    if coll_idx >= len(ax.collections):
        break
    subset_nes = gsealist[gsealist['condition'] == cond]['NES'].values
    ax.collections[coll_idx].set_facecolor(cmap_obj(norm(subset_nes)))
ax.yaxis.tick_right()
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gradient), ax=ax, location = 'left')
cbar.set_label('Enrichment')
plt.title("GSEA Cell Compartment, Top " + str(top) + ' by FDR', fontsize = '30')
ax.set_xticklabels([label_map[c] for c in cat_order])
ax.tick_params(axis='x', labelsize=20, rotation=90)
ax.tick_params(axis='y', labelsize = 30)
ax.set_xlabel('')
ax.set_ylabel('')
ax.get_legend().remove()
# plt.show()
plt.close()
gseaplot.savefig(os.path.join(fig_out, 'protacs_10_de_gsea_chemicalseries_1and10uM.pdf'))

# %% LFC line-plots
matrix_0p1 = pd.read_csv(os.path.join(data_dir, 'Cluster_LFCxPval_0p1uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T
matrix_1 = pd.read_csv(os.path.join(data_dir, 'Cluster_LFCxPval_1uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T
matrix_10 = pd.read_csv(os.path.join(data_dir, 'Cluster_LFCxPval_10uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T
title = 'GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX'
toplot = gsea_df.loc[title]['Lead_genes'][1].split(';')
toplot = gsea_df.loc[title]['Lead_genes'][1].split(';')

# Extract contrast data across all tested concentrations to build dose-response curve
def format_df(df):
    """
    Pick target rows from a concentration matrix and remap indexes. 
    """    
    df1 = df.iloc[[5,14, 11, 0]]
    df1.index = [14, 9, 6, 1]
    return df1

matrix_0p1 = format_df(matrix_0p1)
matrix_0p1['conc'] = 0.1
matrix_1 = format_df(matrix_1)
matrix_1['conc'] = 1
matrix_10 = format_df(matrix_10)
matrix_10['conc'] = 10
m = [matrix_0p1, matrix_1, matrix_10]
m = pd.concat(m)
t_matrix = m[toplot]
t_matrix['conc'] = m['conc']

# Plot dose-response lineplots 
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'
lineplot2 = plt.figure(figsize = (8,10))
plt.errorbar(t_matrix.loc[14]['conc'], t_matrix.loc[14].drop('conc', axis = 1).mean(axis = 1), t_matrix.loc[14].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'orangered', label = 'CRBN-lenalinomide-5N (cluster 14)', linewidth = 2.5)
plt.errorbar(t_matrix.loc[9]['conc'], t_matrix.loc[9].drop('conc', axis = 1).mean(axis =1), t_matrix.loc[9].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'lightgreen', label = 'CRBN-thalidomide-6N (cluster 9)', linewidth = 2.5)
plt.errorbar(t_matrix.loc[6]['conc'], t_matrix.loc[6].drop('conc', axis = 1).mean(axis =1), t_matrix.loc[6].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'aqua', label = 'CRBN-thalidomide-5N (cluster 6)', linewidth = 2.5)
plt.errorbar(t_matrix.loc[1]['conc'], t_matrix.loc[1].drop('conc', axis = 1).mean(axis = 1), t_matrix.loc[1].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'darkblue', label = 'VHL-amide (cluster 1)', linewidth = 2.5)
# plt.legend(fontsize = 26)
plt.ylabel('Averaged DE (' + str(t_matrix.shape[1] )+ ' Proteins)')
plt.xlabel('Concentration \u03bcM')
plt.title(title)
# plt.show()
plt.close()
lineplot2.savefig(os.path.join(fig_out,'protacs_10_' + title.lower() + '_lineplot.pdf'))
