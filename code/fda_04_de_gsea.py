#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: fda_04_de_gsea.py
Description: Perform preranked GSEA on FDA test-set limma contrasts vs DMSO, 
             generate enrichment strip plots, and export leading-edge protein heatmaps.

Author: Shaon Basu
Date: 2025-09-16

Inputs
------
- data/FDA_LimmaMatrix_250304a.csv
- data/FDA_adjLimmaMatrix_250304a.csv
- data/c5.go.v2023.2.Hs.symbols.gmt

Outputs
-------
- figures/fda_04_gsea_targets.pdf
- figures/fda_04_leadingedge_oxdetox_MTX.pdf

Requirements
------------
Python >= 3.8  
Dependencies: pandas, numpy, matplotlib, seaborn, gseapy

"""

# %% Import packages
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
import gseapy as gp
import os

# Set relative paths
data_out = os.path.dirname(__file__)
figure_out = os.path.join(data_out, '..', 'figures')
data_out = os.path.join(data_out, '..', 'data')
workflow = 'fda_04'

# %% Read in limma contrast data vs DMSO from FDA test set & clean
t_matrix = pd.read_csv(os.path.join(data_out,'FDA_LimmaMatrix_250304a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T
t_matrix.index = t_matrix.index.str.replace('Drug_', '')
t_matrix.index = t_matrix.index.str.replace(' - DMSO', '')

# %% GSEA on limma contrasts
def run_gsea(x):
  """
  Run a preranked GSEA analysis.

  Parameters
  ----------
  x : pandas.DataFrame or pandas.Series
      Ranked gene list (e.g., gene symbols with scores).

  Returns
  -------
  gsea_results : gseapy.gsea.GSEA
      GSEA results object with enrichment statistics.
  """
  gsea_results = gp.prerank(rnk=x, 
                          gene_sets=os.path.join(data_out, 'c5.go.v2023.2.Hs.symbols.gmt'), 
                          outdir=data_out,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)
  gsea_df = gsea_results.res2d
  return gsea_df

# Perform GSEA on limma differential expression contrasts vs DMSO
df1 = run_gsea(t_matrix.loc['Clotrimazole'])
df2 = run_gsea(t_matrix.loc['Doxorubicin..Adriamycin..HCl'])
df3 = run_gsea(t_matrix.loc['Epirubicin.HCl'])
df4 = run_gsea(t_matrix.loc['Methotrexate'])
df7 = run_gsea(t_matrix.loc['Ibuprofen.'])
df9 = run_gsea(t_matrix.loc['Doripenem.Hydrate'])

# %% Plot GSEA by NES and FDR cutoffs
gsea_dict = {
    "Clotrimazole": df1,
    "Doxorubicin": df2,
    "Epirubicin": df3,
    "Methotrexate": df4,
    "Ibuprofen": df7,
    "Doripenem": df9}
df_store = pd.DataFrame()  
cat_order = ['Ibuprofen', 'Methotrexate','Epirubicin', 'Doxorubicin', 'Doripenem', 'Clotrimazole']
filter1 = 0.01 #FDR cutoff
filter2 =  0 #NES cutoff
top_pathways = 3

# Create top GSEA table
for idx, df in gsea_dict.items():
    df_subset = df.loc[df['FDR q-val'] < filter1]
    df_subset = df_subset.loc[abs(df_subset['NES']) > filter2]
    df_subset.index = df_subset['Term']
    df_subset['condition'] = idx
    forward = df_subset.sort_values(by = 'NES')
    forward = forward.head(top_pathways)
    reverse = df_subset.sort_values(by = 'NES', ascending = False)
    reverse = reverse.head(top_pathways)
    combined = pd.concat([reverse, forward])
    df_subset = combined
    df_store = pd.concat([df_store,df_subset], axis =0)

# Rename pathways in 'df_store' to be readable
# Split and take substring
df_store['Term'] = [idx.split('_', 1)[1] if '_' in idx else idx for idx in df_store['Term']]
# Term filter
df_store['Term'] = ['_'.join(idx.split('_')[-9:]) if idx.count('_') >= 6 else idx for idx in df_store['Term']]
# Convert to lowercase 
df_store['Term'] = df_store['Term'].str.lower().str.replace('_', ' ')

# Setup GSEA result plot
norm = TwoSlopeNorm(vmin=df_store['NES'].min(), vcenter=0, vmax=df_store['NES'].max())
gradient = 'RdBu_r'
custom_palette = ['blue','red'] 
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'
# Create a stripplot (Pythonic dotplot)
df_store = df_store.rename_axis('C5')
gseaplot = plt.figure(figsize=(5, 20))
sns.stripplot(x="condition", y="Term", data=df_store, hue="NES",  order = cat_order,
              palette=gradient, jitter=False, size=30, dodge=False, legend = False)
ax = plt.gca()
ax.yaxis.tick_right()
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gradient), ax=ax, location = 'top')
cbar.set_label('Enrichment')
ax.set_ylabel("Gene Set Enrichment Analysis \nTop " + str(top_pathways) + " up down, FDR < " + str(filter1), fontsize = '30')
ax.tick_params(axis='x', labelsize=30, rotation=90)  
ax.tick_params(axis='y', labelsize = 30)  
ax.set_xlabel('') 
#plt.show()
gseaplot.savefig(os.path.join(figure_out, f'{workflow}_gsea_targets.pdf')) # Export dotplot plot of GSEA on limma contrasts in FDA-test set

# %% Export GSEA protein matrix within leading edge for a target pathway as a heatmap
adj_matrix = pd.read_csv(os.path.join(data_out,'FDA_adjLimmaMatrix_250304a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T
adj_matrix.index = adj_matrix.index.str.replace('Drug_', '')
adj_matrix.index = adj_matrix.index.str.replace(' - DMSO', '')

mtx = df4.copy()
term = 'GOBP_CELLULAR_OXIDANT_DETOXIFICATION' # Selected term (target pathway
mtx = mtx[mtx['Term']==term]['Lead_genes'].iloc[0].split(';')
qval = str(round(df4[df4['Term']==term]['FDR q-val'].iloc[0],3))
f_matrix = adj_matrix[mtx]
tolabel = ['Methotrexate', 'Doxorubicin..Adriamycin..HCl', 'Pirarubicin', 'Epirubicin.HCl', 'Fulvestrant','Clotrimazole']
f_matrix = f_matrix.loc[tolabel]
f_matrix = f_matrix.T
f_matrix.columns = ['Methotrexate', 'Doxorubicin', 'Pirarubicin', 'Epirubicin', 'Fulvestrant','Clotrimazole']
sns.set(font_scale = 1)
g = sns.clustermap(f_matrix,figsize=(10,10),cmap='bwr', row_cluster = False, col_cluster=False, cbar_pos=(0.1, .2, .03, .1))
g.ax_heatmap.set_title(term.replace('_',' ') + '\nLeading Edge Subset in MTX (FDR < ' + qval + ')\n', fontsize = 20)
g.ax_heatmap.set_ylabel('Proteins', fontsize = 16, rotation = 90)
g.cax.set_title('Moderated t-test \n-log10(adjpval) \nx sign(LFC)')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize = 16)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize = 16)
plt.savefig(os.path.join(figure_out, f'{workflow}_leadingedge_oxdetox_MTX.pdf')) # Export leading edge member protein heatmap (target GO term) for Methotrexate contrast vs DMSO
plt.close()
#plt.show()
