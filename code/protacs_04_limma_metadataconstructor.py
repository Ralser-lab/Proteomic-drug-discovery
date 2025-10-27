#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_04_limma_metadataconstructor.py
Description: Format HBD screen proteome matrix and metadata for compatibility 
             with the R limma package.

Author: Shaon Basu
Date: 2025-09-16

Inputs
------
- data/SB_PROTACs_prmatrix_filtered_5_imputed_50_ltrfm_batched_summarized_240314.tsv
- data/SB_PROTACs_metadata_240611a.tsv

Outputs
-------
- data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
- data/SB_PROTACs_metadata_clustered_240611a.tsv

Requirements
------------
Python >= 3.8
Dependencies: pandas

"""
# %% Import packages
import pandas as pd
import os

# Set relative paths
dir = os.path.dirname(__file__)
data = os.path.join(dir, '..', 'data')
filepath2 = os.path.join(data, 'SB_PROTACs_metadata_240611a.tsv')
export_path = data
import_path = data

# %% Load summarized proteome and metadata from HBD screen
pasef_summarized = pd.read_csv(os.path.join(import_path, 'SB_SpeedyPasef_prmatrix_plateswap_filtered_5_imputed_50_ltrfm_batched_summarized_251027.tsv'),
                     decimal=',', 
                     delimiter=';', 
                    index_col = 0)

pasef_metadata = pd.read_csv(filepath2, 
                         index_col = 0)

# %% Extract metadata columns of interest for limma 
pasef_metadata_clustered = pasef_metadata.drop(pasef_metadata.columns[:28], axis = 1)

# %% Export for limma
pasef_summarized.T.to_csv(os.path.join(export_path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                        header = True, index=True)

pasef_metadata_clustered.to_csv(os.path.join(export_path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv'),
                            header = True, index = True)

