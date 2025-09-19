#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_11_de_stringnetworkenrich_250305a.py
Description:
    Process STRING network enrichment results (GO pathways) for chemical 
    series 4 and 14, identify the top enriched pathways by false discovery 
    rate (FDR), and visualize normalized enrichment scores (NES) using a 
    stripplot. The script filters enriched terms by an FDR cutoff and 
    highlights condition-specific pathway enrichment.

Author: Shaon Basu
Date: 2025-09-19

Inputs
------
- data/Cluster4_enrichment.Component.tsv
- data/Cluster14_enrichment.Component.tsv

Outputs
-------
- figures/protacs_11_STRINGenrich_PROTACs_OnandOfftarg.pdf

Requirements
------------
Python >= 3.8
Dependencies: pandas, numpy, matplotlib, seaborn, os

"""

# %% Import modules
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
import numpy as np
import os

# Set relative paths
dir_main = os.path.dirname(__file__)
filepath2 = os.path.join(dir_main, '..', 'data/')
filepath3 = os.path.join(dir_main, '..', 'figures')

# Load string network enrichment values (FDR, NES)
# https://string-db.org/, queried 10 uM DE profiles from 
# chemical series 4 & 14 (vs DMSO [limma]) for GO pathways
TXNreact = pd.read_csv(filepath2 + 'Cluster4_enrichment.Component.tsv', sep = '\t')
L5Nreact = pd.read_csv(filepath2 + 'Cluster14_enrichment.Component.tsv', sep = '\t')

# %% Create enrichment stripplot (dotplot of NES values)
nx_dict = {
    "Nuclear PROTAC": TXNreact,
    "Lenalinomide 5N": L5Nreact,
}
cat_order = ['Nuclear PROTAC', 'Lenalinomide 5N']

def top5(df):
    """
    Return the top 5 rows with the lowest false discovery rate.
    """
    df = df.sort_values('false discovery rate', ascending = True)
    df = df.iloc[0:5]
    return df

for n, df in nx_dict.items():
    df = top5(df)
    nx_dict[n] = df

# Create top GSEA table
dogma = pd.DataFrame()
filter = 0.0001 #FDR cutoff
for idx, df in nx_dict.items():
    df_subset = df.loc[df['false discovery rate'] < filter]
    df_subset.index = df_subset['term description']
    if df_subset['direction'][0] == 'top':
        df_subset['enrichment score'] = -(df_subset['enrichment score'])
    df_subset['condition'] = idx
    dogma = pd.concat([dogma,df_subset], axis =0)

# Rename pathways to be readable
dogma.loc[dogma['enrichment score'] == dogma['enrichment score'].min(), 'enrichment score'] = dogma['enrichment score'].sort_values()[1]
norm = plt.Normalize(vmin=dogma['enrichment score'].min(), vmax=dogma['enrichment score'].max())
gradient = 'RdBu_r'
custom_palette = ['blue','red'] 
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'

# Create a stripplot (dot plot)
gseaplot = plt.figure(figsize=(10, 10))
sns.stripplot(x="condition", y="term description", data=dogma, hue="enrichment score", order = cat_order,
              palette=gradient, jitter=False, size=30, dodge=False, legend = False)
ax = plt.gca()
ax.yaxis.tick_left()
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gradient), ax=ax)
cbar.set_label('Enrichment')
plt.title("PPI enrichment, \nTop 5 Pathways (FDR, NES) " + str(filter), fontsize = '30')
ax.tick_params(axis='x', labelsize=18)  
ax.tick_params(axis='y', labelsize = 18)  
ax.set_xlabel('')  
ax.set_ylabel('') 
#plt.show()
gseaplot.savefig(os.path.join(filepath3, 'protacs_11_STRINGenrich_PROTACs_OnandOfftarg.pdf'))
