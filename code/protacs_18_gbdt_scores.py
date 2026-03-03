#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_18_gbdt_scores.py
Description:
    Post-hoc analysis of xgb model outputs (toxicity scores) with a focus on 
    series 15 compounds and the lenalidomide 5N series analogues. The script 
    loads exported model predictions (Rplot_Figure4.csv) and trained XGBoost 
    models, and generates supplementary tables and figures:

Author: Shaon Basu
Date: 2025-09-30

Inputs
------
- data/Rplot_Figure4.csv
- scoring_models/protacs_16_xgb_second-pass-model.json
- scoring_models/protacs_16_xgb_second-pass-calibrated-model.pkl.

Outputs
-------
- data/NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv
- data/analog_dataout.csv
- figures/protacs_17_series15_toxscores.pdf
- figures/protacs_17_analogues_signature_PCA.pdf
- figures/protacs_17_analogues_signature_barplot.pdf

Requirements
------------
Python >= 3.8
Dependencies:
    - pandas, numpy, matplotlib, seaborn, scikit-learn, xgboost, joblib

"""
# %% Import packages
import pandas as pd
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from device_gradientboostingmachine import GBDT

# Set relative paths
GBDT.configure_font()
workflow = 'protacs_18'
path = os.path.dirname(__file__)
inputout = os.path.join(path, '..', 'data')
figout = os.path.join(path, '..','figures')
modelvalues_w = pd.read_csv(os.path.join(inputout,'Rplot_Figure4_protacs_16.csv'), index_col = 0)
modelvalues_n = pd.read_csv(os.path.join(inputout,'Rplot_Figure4_protacs_17.csv'), index_col = 0)
cols_only_in_w = modelvalues_w.columns.difference(modelvalues_n.columns)
modelvalues = modelvalues_n.join(modelvalues_w[cols_only_in_w]).iloc[:,5:]
cols_df = modelvalues.columns.difference(modelvalues_n)
modelvalues_f = pd.merge(modelvalues_n.iloc[:,:5], modelvalues, left_index=True, right_index=True)
modelvalues_f.to_csv(os.path.join(inputout, 'Rplot_Figure4.csv'))
softvote_df = pd.DataFrame()
softvote_df = softvote_df.assign(
    probability = (modelvalues_w.Probability + modelvalues_n.Probability)/2,
    drug_id = modelvalues_n.Drug,
    cluster = modelvalues_n.Cluster)

# %%                 

# Export Table S3
tableS3 = softvote_df.copy().rename(columns = {
    'probability': 'Toxic Probability',
    'drug_id': 'Drug',
    'cluster': 'Cluster'
})


# Map internal AZ IDs to friendly series names
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14183816').values, 'Drug'] = 'Compound 1'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14196658').values, 'Drug'] = 'Compound 2'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14197166').values, 'Drug'] = 'Compound 3'

# Quick  check of relabeling (no assignment; used when debugging)
tableS3.loc[tableS3['Drug'].isin(['Compound 1', 'Compound 2', 'Compound 3'])]

# Order by cluster and export
tableS3.sort_values('Cluster', inplace=True)
tableS3

# Add inptu weights to probabilities
tableS3 = pd.merge(tableS3, modelvalues, left_index=True, right_index=True)
tableS3.to_csv(os.path.join(inputout, 'NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv'))

# %% Cluster 15 probability plot
# Cluster 15 proba
cluster15_df = softvote_df.loc[softvote_df.cluster==15].drop('drug_id', axis = 1)
# Other chemotype controls proba
control_ids = ['Moxifloxacin', 'Troglitazone', 'Ibuprofen', 'Ambrisentan', 'Ticrynafen']
control_df = softvote_df.loc[softvote_df.drug_id.isin(control_ids)].reset_index().set_index('drug_id').drop('index', axis = 1)
# Get cluster 15 and control chemotypes in one filter
combined_df = pd.concat([cluster15_df, control_df])
# Merge replicates by core ID (strip trailing _SN) and average probability
values = (
    combined_df.probability
    .groupby(combined_df.index.str.rsplit('_', n=1).str[0])
    .mean()
    .sort_index()
    .sort_values(ascending=False)
)
# Plot toxicity probabilities (Series 15 vs. controls)
plt.rcParams['axes.labelsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
plt.rcParams['xtick.labelsize'] = '20'
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
plt.figure(figsize=(10, 6))
values.plot(kind='bar')
plt.savefig(os.path.join(figout, f'{workflow}_series15_toxscores.pdf'))
# plt.show()
plt.close()
values = list(values.index.values)  # keep original indices if needed downstream

# %% AZ14183816 chemical series (analogues) scores
# Select specific analogues; align to modelvalues and append probability
cpd_targets = 'AZ14183816-005, AZ14183816-006, AZ14196658-003, AZ14197166-003, AZ14197166-004, AZ14197166-005'
cpd_targets = cpd_targets.replace('-', '_').replace(' ', '').split(',')
cpd_targets = modelvalues_f.iloc[:, 5:].assign(probability=modelvalues_f['Probability']).loc[cpd_targets]
cpd_targets['analogue'] = cpd_targets.index.str.split('_').str[0]

# %% PCA on analogues
analog = cpd_targets.groupby('analogue').mean()          # collapse replicates to analogue means
analog.copy().drop('probability', axis=1, inplace=True)  # ensure side-effect free in-place drop
analog.to_csv(os.path.join(inputout, 'analog_dataout.csv'))

def norm_and_plot(dataframe):
    """
    Z-score columns for PCA; quick internal barplot preview
    helper function (optional).

    """
    z = (dataframe - dataframe.mean()) / dataframe.std()
    # z.plot(kind='bar')
    return z

# Exclude probability column from PCA
df = norm_and_plot(analog.loc[:, analog.columns != 'probability'])

# Fit PCA on z-scored features
pca = PCA()
pca.fit(df)

# Build loadings and component scores tables
loadings = pd.DataFrame(
    pca.components_.T,
    columns=[f'PC_{i+1}' for i in range(len(df.index))],
    index=df.columns
)
components = pd.DataFrame(
    pca.fit_transform(df),
    columns=[f'PC_{i+1}' for i in range(len(df.index))],
    index=df.index
)

# Plot PCA scatter with biplot arrows
plt.rcParams['axes.labelsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
plt.rcParams['xtick.labelsize'] = '20'
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
ev = pca.explained_variance_ratio_
sns.scatterplot(x='PC_1', y='PC_2', data=components, hue=df.index)
plt.xlabel(f"PC_1 ({ev[0]*100:.2f}%)")
plt.ylabel(f"PC_2 ({ev[1]*100:.2f}%)")

# Scale loadings for arrows & annotate with feature names
loadings_ratio = pca.components_.T * np.sqrt(pca.explained_variance_)
for i, feature in enumerate(loadings.index):
    plt.arrow(0, 0, loadings_ratio[i, 0], loadings_ratio[i, 1],
              color='red', alpha=0.7, head_width=0.05)
    plt.text(loadings_ratio[i, 0] * 2, loadings_ratio[i, 1] * 1.3,
             feature, color='black', ha='center', va='center')
plt.legend(title='Analogue', fontsize='small', title_fontsize='medium')
plt.savefig(os.path.join(figout, f'{workflow}_analogues_signature_PCA.pdf'))
# plt.show()
plt.close()

# %% PC2 top features
# Pick features with strongest PC2 loadings, then plot analogue means
pca_features = loadings.sort_values(by='PC_2').index[0:4]
etc_features = ['NDUFA5', 'NDUFA4', 'CYC1']  # override with Complex I panel

barplot_df = analog[etc_features]

#barplot_df = (barplot_df - barplot_df.mean())/barplot_df.std()
barplot_df.plot(kind='bar')
plt.savefig(os.path.join(figout, f'{workflow}_analogues_signature_barplot.pdf'))
plt.ylim([-0.4,1.2])
# plt.show()
plt.close()
# %%
