# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import sys
sys.path.append('/Users/shaon/Desktop/PROTACS/50-0121_PROTAC_HCC_Lib/Code/')
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import preprocessingdevice
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import leaves_list
import os

# Formats metadata and protein expresison matrix for compatibility with limma package (R)

# %% Set rel working directories

dir = os.path.dirname(__file__)

data = os.path.join(dir, '..', 'data')

# %% Load Datasets

filepath2 = os.path.join(data, 'SB_PROTACs_metadata_240611a.tsv')

export_path = data

import_path = data

pasef_summarized = pd.read_csv(os.path.join(import_path, 'SB_PROTACs_prmatrix_filtered_5_imputed_50_ltrfm_batched_summarized_240314.tsv'),
                     decimal=',', 
                     delimiter=';', 
                    index_col = 0)

pasef_metadata = pd.read_csv(filepath2, 
                         index_col = 0)

# %% Extract columns of interest for limma

pasef_metadata_clustered = pasef_metadata.drop(pasef_metadata.columns[:28], axis = 1)

# %% export for Limma

pasef_summarized.T.to_csv(os.path.join(export_path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                        header = True, index=True)

pasef_metadata_clustered.to_csv(os.path.join(export_path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_metaforlimma_240611a.tsv'),
                            header = True, index = True)

