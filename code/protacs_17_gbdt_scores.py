#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_23_gbdt_scores.py
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
import joblib
import os
import xgboost as xgb
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set relative paths
path = os.path.dirname(__file__)
inputout = os.path.join(path, '..', 'data')
modelout = os.path.join(path, '..', 'scoring_models')
figout = os.path.join(path, '..','figures')
modelvalues = pd.read_csv(os.path.join(inputout,'Rplot_Figure4.csv'), index_col = 0)

# Load trained models
finalmodel = xgb.XGBClassifier()
finalmodel.load_model(os.path.join(modelout, 'protacs_16_xgb_second-pass-model.json'))
calibratedmodel = joblib.load(os.path.join(modelout, 'protacs_16_xgb_second-pass-calibrated-model.pkl'))

# %% Export Table S3
tableS3 = modelvalues.copy()
tableS3.drop(labels=['Actual', 'Predicted'], axis=1, inplace=True) 
tableS3.columns.values[0] = 'Toxic Probability' 

# Map internal AZ IDs to friendly series names
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14183816').values, 'Drug'] = 'Compound 1'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14196658').values, 'Drug'] = 'Compound 2'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14197166').values, 'Drug'] = 'Compound 3'

# Quick  check of relabeling (no assignment; used when debugging)
tableS3.loc[tableS3['Drug'].isin(['Compound 1', 'Compound 2', 'Compound 3'])]

# Order by cluster and export
tableS3.sort_values('Cluster', inplace=True)
tableS3.to_csv(os.path.join(inputout, 'NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv'))

# %% Series 15 and control predictions
cluster15 = modelvalues.loc[modelvalues['Cluster'] == 15].sort_values('Probability')
targets = ['Moxifloxacin', 'Troglitazone', 'Ibuprofen', 'Ambrisentan', 'Ticrynafen']
controls = modelvalues.loc[modelvalues['Drug'].isin(targets)].sort_values('Probability')
controls.index = controls['Drug']  # label with drug names for plotting
toplot = pd.concat([cluster15, controls])

# Merge replicates by core ID (strip trailing _SN) and average probability
values = (
    toplot['Probability']
    .groupby(toplot.index.str.rsplit('_', n=1).str[0])
    .mean()
    .sort_index()
    .sort_values()
)

# Plot toxicity probabilities (Series 15 vs. controls)
plt.rcParams['axes.labelsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
plt.rcParams['xtick.labelsize'] = '20'
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
plt.figure(figsize=(10, 6))
values.plot(kind='bar')
plt.savefig(os.path.join(figout, 'protacs_17_series15_toxscores.pdf'))
# plt.show()
plt.close()

values = list(toplot.index.values)  # keep original indices if needed downstream

# %% AZ14183816 chemical series (analogues) scores
# Select specific analogues; align to modelvalues and append probability
cpd_targets = 'AZ14183816-005, AZ14183816-006, AZ14196658-003, AZ14197166-003, AZ14197166-004, AZ14197166-005'
cpd_targets = cpd_targets.replace('-', '_').replace(' ', '').split(',')
cpd_targets = modelvalues.iloc[:, 5:].assign(probability=modelvalues['Probability']).loc[cpd_targets]
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
sns.scatterplot(x='PC_1', y='PC_2', data=components, hue=df.index)

# Scale loadings for arrows & annotate with feature names
loadings_ratio = pca.components_.T * np.sqrt(pca.explained_variance_)
for i, feature in enumerate(loadings.index):
    plt.arrow(0, 0, loadings_ratio[i, 0], loadings_ratio[i, 1],
              color='red', alpha=0.7, head_width=0.05)
    plt.text(loadings_ratio[i, 0] * 2, loadings_ratio[i, 1] * 1.3,
             feature, color='black', ha='center', va='center')
plt.legend(title='Analogue', fontsize='small', title_fontsize='medium')
plt.savefig(os.path.join(figout, 'protacs_17_analogues_signature_PCA.pdf'))
# plt.show()
plt.close()

# %% PC2 top features
# Pick features with strongest PC2 loadings, then plot analogue means
features = loadings.sort_values(by='PC_2').index[0:4]
features = ['NDUFA5', 'CYC1', 'NDUFA4', 'NDUFB10', 'NDUFA13']  # override with Complex I panel
analog[features].plot(kind='bar')
plt.savefig(os.path.join(figout, 'protacs_17_analogues_signature_barplot.pdf'))
# plt.show()
plt.close()