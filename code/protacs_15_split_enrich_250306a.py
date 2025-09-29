#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_15_split_enrich_250306a.py
Description:
    Process differential expression (DE) results and Reactome enrichment
    analyses for HBD chemical series from the split dataset. The script 
    generates signed ranked gene lists from limma DE outputs on the split, 
    filters pathway enrichment results by false discovery rate (FDR), and
    assembles a summary of top enriched pathways across chemical clusters.

Author: Shaon Basu
Date: 2025-09-19

Inputs
------
- data/Cluster_VHL_Limma_250306a.csv
- data/Cluster_T6N_Limma_250306a.csv
- data/Cluster_T5N_Limma_250306a.csv
- data/Cluster_L5N_Limma_250306a.csv
- data/Cluster_URA_Limma_250306a.csv
- data/Cluster_TXN_Limma_250306a.csv
- data/DEsplit_*_reactome.tsv (Reactome enrichment results per cluster)

Outputs
-------
- data/DE_VHL.csv, DE_T6N.csv, DE_T5N.csv, DE_L5N.csv, DE_URA.csv, DE_TXN.csv
  (signed ranked gene scores per HBD chemical series [recruiter type] vs DMSO)
- data/top5_FDR_reactome.csv
  (summary of enriched Reactome pathways at FDR < 0.01 across series)

Requirements
------------
Python >= 3.8
Packages: pandas, numpy, os

"""

# %% Import modules
import pandas as pd
import numpy as np
import numpy as np
import os

# Set relative paths
base = os.path.dirname(__file__)
outpath = os.path.join(base, '..', 'data')
figures = os.path.join(base, '..', 'figures')

# Load differential expression profiles (chemical series vs DMSO) from split
DE_VHL = pd.read_csv(outpath + '/Cluster_VHL_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)
DE_T6N = pd.read_csv(outpath + '/Cluster_T6N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)
DE_T5N = pd.read_csv(outpath + '/Cluster_T5N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)
DE_L5N = pd.read_csv(outpath + '/Cluster_L5N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)
DE_URA = pd.read_csv(outpath + '/Cluster_URA_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)
DE_TXN = pd.read_csv(outpath + '/Cluster_TXN_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

# Load string network enrichment values (FDR, NES)
# https://string-db.org/, queried 10 uM DE profiles from 
# chemical series in split (listed below) for GO pathways
DE_list = {
    "DE_VHL": DE_VHL,
    "DE_T6N": DE_T6N,
    "DE_T5N": DE_T5N,
    "DE_URA": DE_URA,
    "DE_L5N": DE_L5N,
    "DE_TXN": DE_TXN,
}

def process_for_network(df):
    """
    Create a signed ranking of genes from logFC and p-values.
    Returns a Series of -log10(p-value) weighted by the sign of logFC,
    sorted in descending order.
    """
    df['sign'] = np.sign(df['logFC']) * -np.log10(df['P.Value'])
    ranked = df['sign'].sort_values(ascending=False)
    return ranked

for name, df in DE_list.items():
    ranked = process_for_network(df)
    path = os.path.join(outpath, name + '.csv')
    ranked.to_csv(path)

T5Nreact = pd.read_csv(outpath + '/DEsplit_T5N_reactome.tsv', sep = '\t')
TXNreact = pd.read_csv(outpath + '/DEsplit_TXN_reactome.tsv', sep = '\t')
T6Nreact = pd.read_csv(outpath + '/DEsplit_T6N_reactome.tsv', sep = '\t')
L5Nreact = pd.read_csv(outpath + '/DEsplit_L5N_reactome.tsv', sep = '\t')
VHLreact = pd.read_csv(outpath + '/DEsplit_VHL_reactome.tsv', sep = '\t')
URAreact = pd.read_csv(outpath + '/DEsplit_URA_reactome.tsv', sep = '\t')

nx_dict = {
    "Thalidomide 5N": T5Nreact,
    "Transcription": TXNreact,
    "Dihydrouracyl": URAreact,
    "Lenalinomide 5N": L5Nreact,
    "Thalidomide 6N": T6Nreact,
    "VHL amide": VHLreact
}

def top5(df):
    """
    Return the top 10 rows with the lowest false discovery rate.
    """
    df = df.sort_values('false discovery rate', ascending = True)
    df = df.iloc[0:10]
    return df

for n, df in nx_dict.items():
    df = top5(df)
    nx_dict[n] = df

cat_order = ['Thalidomide 5N','Lenalinomide 5N','VHL amide','Thalidomide 6N', 'Dihydrouracyl','Transcription']
#cat_order = ['Thalidomide 5N','Lenalinomide 5N','VHL amide','Thalidomide 6N','Transcription']
cat_order.reverse() 
dogma = pd.DataFrame()
filter = 0.01 #FDR cutoff

# Assemble a cross-condition summary table of enriched pathways in split and export 
for idx, df in nx_dict.items():
    df_subset = df.loc[df['false discovery rate'] < filter]
    df_subset.index = df_subset['term description']
    if df_subset['direction'][0] == 'top':
        df_subset['enrichment score'] = -(df_subset['enrichment score'])
    df_subset['condition'] = idx
    dogma = pd.concat([dogma,df_subset], axis =0)
dogma.sort_values('enrichment score', ascending = True)
dogma['term description'] = pd.Categorical(dogma['term description'], ordered = True)
dogma.to_csv(outpath + '/top5_FDR_reactome.csv')
