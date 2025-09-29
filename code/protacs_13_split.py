#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_13_split.py
Description:
    Prepare training and test splits for gradient boosting machine learning 
    models of mitochondrial toxicity. The script cleans drug identifiers, 
    formats proteomic log fold change (LFC) data, aligns with IC50 values, 
    and generates a binary toxicity outcome variable. The dataset is then 
    split into training and test sets, and compound metadata is updated 
    with split labels for downstream modeling.

Author: Shaon Basu
Date: 2025-09-19

Inputs
------
- data/Drug_LFCxPval_250305a.csv
- data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
- data/AZcompound_metadata_clustered_240611a.tsv

Outputs
-------
- data/AZcompound_metadata_clustered_split_240611a.tsv

Requirements
------------
Python >= 3.8  
Dependencies: pandas, scikit-learn, os, sys  

"""
# %% Import modules
import pandas as pd
import os
import sys
from sklearn.model_selection import train_test_split

# Set relative paths
path = os.path.dirname(__file__)
path = os.path.join(path, '..', 'data')
out = os.path.join(path, '..', 'figures')

# Cleaning function for data loading
def clean_drug_index(df):
    """
    Clean the DataFrame index by extracting drug IDs and removing extra text.
    """
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

# Load datasets
LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 'Drug_LFCxPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)
expression_matrix =  pd.read_csv(os.path.join(path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T

# DIA-MS input for Machine Learning
matrix = LFC_matrix.copy()

# %% Format predictor and response vars
AZmeta_input= pd.read_csv(os.path.join(path, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)
AZmeta = AZmeta_input.copy()
AZmeta.index = AZmeta.index.str.replace('-','_') # Format index
AZmeta['Gal'] = AZmeta['Gal'].str.replace('>','').astype(float) # Convert IC50 to float
IC50s = pd.DataFrame({'Gal_IC50': AZmeta['Gal']}, dtype = float) # Extract IC50s
IC50s = IC50s.loc[matrix.index] # Match to proteome index
IC50snoNA = IC50s[IC50s['Gal_IC50'].isna()==False] # remove NAs
BinaryTox = pd.DataFrame({ # Make outcome var categorical
    'IC50' : (IC50snoNA['Gal_IC50']<10).astype(int)  
    })

# %% Produce train-test split
X = matrix.loc[BinaryTox.index].copy() 
y = BinaryTox.copy()
Xb_train, Xb_test, yb_train, yb_test = train_test_split(X, BinaryTox, test_size=0.2, random_state=42) # Baseline Split
AZmeta_input.loc[Xb_train.index.str.replace('_','-'), 'GBDT_split_idx'] = 1
AZmeta_input.to_csv(os.path.join(path, 'AZcompound_metadata_clustered_split_240611a.tsv'))
