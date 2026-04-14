
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_09_de_pca.py
Description: Perform PCA and probability-based analyses on differential
             expression profiles from PROTAC and FDA drug screening datasets. 
             The script generates 3D PCA plots, probability mass/density plots, 
             violin/boxplots, runs z-tests, and exports processed summary tables.

Author: Shaon Basu
Date: 2025-09-18

Inputs
------
- data/Cluster_LFCxPval_10uM_250305a.csv
- data/Drug_LFCxadjPval_250305a.csv
- data/AZcompound_metadata_clustered_240611a.tsv
- data/FDA_proba_250304.tsv

Outputs
-------
- figures/protacs_09_de_pca_chemicalseries.pdf
- figures/protacs_09_de_kde_fda_v_hbd_zcentered.pdf
- figures/protacs_09_de_kde_fda_v_hbd_boxplot.pdf
- data/NCB_ProteomeGuidedDiscovery_TableS2_250606a.csv

Requirements
------------
Python >= 3.8
Dependencies: pandas, numpy, seaborn, matplotlib, scikit-learn
Custom: device_summarystatistics

""" 
# %% import modules
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import numpy as np

# Import custom statistics module
import sys
sys.path.append(os.path.join(os.path.dirname(__file__)))
import device_summarystatistics

# Set relative paths
dir1 = os.path.dirname(__file__)
filepath2 = os.path.join(dir1, '..', 'data')
filepath3 = os.path.join(dir1, '..', 'figures')
workflow = 'protacs_09'

# Load differential expression matrix for HBD library screen (chemical series vs DMSO)
t_matrix = pd.read_csv(filepath2 + '/Cluster_LFCxPval_10uM_250305a.csv',
                       delimiter=';', decimal=',', index_col=0, header=0).T

# %% Perform principal component analysis on differential expression profiles (PCA on chemical series vs DMSO)
pca = PCA(n_components=3)
LFC_pca = pca.fit_transform(t_matrix)

# Create a DataFrame with PCA results
df_pca = pd.DataFrame({'PC1': LFC_pca[:, 0], 'PC2': LFC_pca[:, 1], 'PC3': LFC_pca[:, 2]},
                      index=t_matrix.index)

# Add labels and colors for plotting
df_pca['labels'] = df_pca.index.str.replace(' - Cluster2_DMSO','')
df_pca['labels'] = df_pca['labels'].str.replace('Cluster2_','')
unique_labels = df_pca['labels'].unique()
explained_variance_ratios = pca.explained_variance_ratio_

#  Calculate PC1 distance
df_pca['labels'] = df_pca['labels'].astype('float64')
df_pca.sort_values(by = ['PC1'], inplace=True)
cmap = plt.get_cmap('jet')
color_dict = {drug: cmap(i / len(df_pca.index)) for i, drug in enumerate(df_pca['labels'])}
df_pca['color'] = df_pca['labels'].map(color_dict)

# Manually remap one-hot encoded labels for dendrograms based on 
# PC1 distance after PCA on differential expression profiles of series vs DMSO
# Retroactive mapping in the following scripts for one-hot encoding & dendrogram plotting:
#  - protacs_02_azmetadata_cluster_250306a.py
#  - protacs_03_azmetadata_dendrogram_250306a.R
df_pca['labels'] = df_pca['labels'].astype('string')
clusterlabs = {'1.0': 'AR-VHL-Other',
               '2.0': 'AR-T6N-Indole',
               '3.0': 'AR-Other-Other',
               '4.0': 'Txn-VHL/6N-Other',
               '5.0': 'AR-VHL-Indole',
               '6.0': 'AR-T5N-Other',
               '7.0': 'Txn-Other-Other',
               '8.0': 'AR-T6N/VHL-Pip',
               '9.0': 'AR-T6N-Other',
               '10.0': 'AR-T5N-Pip',
               '11.0': 'AR-DHU-Pip',
               '12.0': 'AR-Other-Pip',
               '13.0': 'Non PROTAC',
               '14.0': 'AR-L5N-Other',
               '15.0': 'AR-L5N-Pip'} 
df_pca['clusterlabs'] = df_pca['labels'].map(clusterlabs)

# Set plot parameters
plt.rcParams['axes.labelsize'] = '30'
plt.rcParams['axes.titlesize'] = '27'
plt.rcParams['xtick.labelsize'] = '26'
plt.rcParams['ytick.labelsize'] = '26'
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot  
scatter = ax.scatter(df_pca['PC1'], df_pca['PC3'], df_pca['PC2'], c=df_pca['color'], s=200, alpha=1)

# Annotate each point with the corresponding target label
for i in range(len(df_pca)):

    ax.text(df_pca['PC1'][i] - 10, df_pca['PC3'][i] + 10, df_pca['PC2'][i] - 3, df_pca['labels'][i], fontsize=15)

# Set labels
ax.set_xlabel(f"PC1 ({explained_variance_ratios[0]*100:.2f}%)")
ax.set_ylabel(f"PC3 ({explained_variance_ratios[1]*100:.2f}%)")
ax.set_zlabel(f"PC2 ({explained_variance_ratios[2]*100:.2f}%)")

# Remove grid lines and background pane
ax.grid(False)
ax.zaxis.pane.fill = False

# Create custom legend
df_pca['labels'] = df_pca['labels'].astype('float64')
df_pca.sort_values(by = ['PC1'], inplace = True)
label_colors = df_pca.drop_duplicates('clusterlabs').set_index('clusterlabs')['color']
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=12) for color in label_colors]
legend = ax.legend(handles, label_colors.index, bbox_to_anchor=(1.05, 1), loc='upper left', title='Drug Target', title_fontsize=26, fontsize=26)
plt.title('Principal Component Analysis')
# plt.show()
# Export PCA on DE profiles in HBD screen (chemical series vs DMSO)
fig.savefig(os.path.join(filepath3,f'{workflow}_de_pca_chemicalseries.pdf'))
plt.close(fig)

# %% Probability mass function (PMF) plots based on DE counts
# Cleans index of differential expression profiles to align with metadatafile
def clean_drug_index(df):
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

# Clean and import data
LFCxadjPval_matrix = clean_drug_index(pd.read_csv(os.path.join(filepath2, 'Drug_LFCxadjPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

# Clean and import metadata
AZmeta= pd.read_csv(os.path.join(filepath2, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)
AZmeta.index = AZmeta.index.str.replace('-','_')

# Condition on adjusted p-values based on alpha cutoff
mask = abs(LFCxadjPval_matrix) > -np.log10(0.05)

# Counts after conditioning on DE profiles
counts = mask.sum(axis = 1).sort_values()/mask.shape[1]
ontargetHBD = counts.loc[AZmeta['Drug_Type'] == 'Txn-PROTAC']
arHBD = counts.loc[AZmeta['Drug_Type'] == 'AR-PROTAC']


# %% Probability density function (PDF) plots based on DE counts
# Load diffexp probabilities from FDA drug screen dataset
counter = pd.read_csv(os.path.join(filepath2, 'FDA_proba_250304.tsv'), index_col = 0).iloc[:,0]

# Combine FDA and HBD drug screening profiles conditioned on DE alpha cut-off
df = pd.DataFrame({
    'FDA approved drugs': counter,
    'AR-directed HBDs': arHBD,
    'on-target HBDs': ontargetHBD
})

# Define custom colors
colors = ['grey', 'red', 'blue']

# Plot density for each column with specified colors
plt.figure(figsize=(8,6))  # Adjust figure size for better visualization
for column, color in zip(df.columns, colors):
    sns.kdeplot(df[column].dropna(), clip=(0, None), label=column, fill=True, alpha=0.5, color=color)

# Add legend and show plot
plt.legend()
plt.title("Proteome response to \ndegrader library")

# Save probability density function as a PDF
plt.savefig(os.path.join(filepath3, f'{workflow}_de_kde_fda_v_hbd_zcentered'))
#plt.show()
plt.close()

# %% Plot probability density functions as boxplot
from matplotlib import pyplot as plt
df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df.items()]))

# Prepare data by dropping NaN values for each column
data_no_nan = [df[col].dropna().values for col in df.columns]

# Create and save the boxplot
plt.figure(figsize=(8, 6))
plt.violinplot(data_no_nan)
plt.boxplot(data_no_nan, labels = df.columns)
plt.title("Proteome response to \ndegrader library")
plt.ylabel("Values")
plt.savefig(os.path.join(filepath3, f'{workflow}_de_kde_fda_v_hbd_boxplot'))
#plt.show()
plt.close()

# Perform z-test with CDF
device_summarystatistics.z_test(df.iloc[:,0], df.iloc[:,1])

# %% Export datasets
# Create df with count info by drug
count_matrix = pd.DataFrame({'deg count': mask.sum(axis = 1),
'total proteins': mask.shape[1], 
'probability': mask.sum(axis=1)/mask.shape[1]
})

# Merge with metadata
tableS2 = pd.merge(AZmeta.loc[count_matrix.index][['Drug ID','Drug_Type', 'Binned_ligase','Binned_Target','Dend','Gal']], 
                count_matrix, 
                left_index = True,
                right_index = True,
                how = 'inner')


# Rename df attributes for interpretability
tableS2.columns = ['Drug Name', 'Drug Type', 'Recruiter', 'Binder', 'Cluster', 'Gal IC50', 'DEG count', 'Total Proteome', 'Probability']
tableS2.replace(['AR-PROTAC','Non-PROTAC','Txn-PROTAC'], ['AR-HBD','not-HBD','ontarget-HBD'], inplace = True)

# Sort by chemical series (remapped on PC1)
tableS2.sort_values(by = 'Cluster', inplace = True)

# Merge with LFC-signed -log10 adjusted p-value matrix
tableS2full = pd.merge(tableS2, LFCxadjPval_matrix, left_index = True, right_index = True)
tableS2.index.name = 'Unique Identifier'

# Retroactive analogue pseudonyms
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14183816').values, 'Drug Name'] = 'Compound 1'
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14196658').values, 'Drug Name'] = 'Compound 2'
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14197166').values, 'Drug Name'] = 'Compound 3'

# Check analogue edit
tableS2.loc[tableS2['Drug Name'].isin(['Compound 1', 'Compound 2', 'Compound 3'])]

# Saveout
tableS2full.to_csv(os.path.join(filepath2, 'NCB_ProteomeGuidedDiscovery_TableS2_250606a.csv'))

