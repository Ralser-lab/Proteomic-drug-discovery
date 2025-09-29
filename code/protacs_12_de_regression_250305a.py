#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_12_de_regression_250305.py
Description:
    Analyze proteomic log fold change (LFC) profiles from the HBD proteome 
    screen in relation to mitochondrial toxicity (IC50, galactose assay). 
    The workflow extracts pathway genes from STRING enrichment data, associates 
    per-series LFCs with IC50 values, and performs univariate linear regression to 
    identify genes predictive of pIC50. Visualization includes IC50 barplots, 
    regression diagnostics, and gene-level observed vs predicted plots.

Author: Shaon Basu
Date: 2025-09-19

Inputs
------
- data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
- data/Cluster_LFC_10uM_250305a.csv
- data/Cluster14_enrichment.Component.tsv
- data/Cluster4_enrichment.Component.tsv
- data/AZcompound_metadata_clustered_240611a.tsv

Outputs
-------
- figures/protacs_12_barplot_IC50.pdf
- figures/protacs_12_<first_pathway> models.pdf
- figures/protacs_12_regression_p_values.pdf
- figures/protacs_12_OLS_<gene>.pdf
- figures/protacs_12_OLS_<gene>_pred.pdf

Requirements
------------
Python >= 3.8  
Dependencies: pandas, numpy, matplotlib, statsmodels, os

"""
# %% Import modules
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# Set relative paths
dir_main = os.path.dirname(__file__)
filepath = os.path.join(dir_main, '..', 'data/')
outpath = os.path.join(dir_main, '..', 'figures/')

# Load expression matrix, differential expression profiles, and enrichment data from HBD proteomic screen. 
expression_matrix =  pd.read_csv(filepath + 'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv',
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T
LFC_matrix = pd.read_csv(filepath + 'Cluster_LFC_10uM_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T
string14 = pd.read_csv(filepath + 'Cluster14_enrichment.Component.tsv', index_col = 1, sep = '\t')
string4 = pd.read_csv(filepath + 'Cluster4_enrichment.Component.tsv', index_col = 1, sep = '\t')

# %% IC50-Proteome regression workflow
# Performs: Pathway gene extraction, IC50 analysis, and regression modelling
pathway1 = list(string14.index[0:5]) + list(string4.index[0:5])  # take first 5 pathway names from string14 and string4

# pathway1 can be manually overridden by setting a single pathway name instead:
# pathway1 = ['Inner mitochondrial membrane protein complex']
# pathway1 = ['Oxidoreductase complex']
# pathway1 = ['Mitochondrial ribosome']
# pathway1 = ['Mitochondrial respirasome']
# pathway1 = ['Respiratory chain complex']
# pathway1 = ['Chromosomal region']
# pathway1 = ['Nuclear chromosome']
# pathway1 = ['MCM complex']
# pathway1 = ['Cytosolic ribosome']
# pathway1 = ['Replication fork']

def extract_genes(df, pathways, matrix, boundary):
    """
    Extract gene list for given pathways, then filter for genes 
    with mean expression above boundary.
    """
    all_genes = set()

    for pathway in pathways:
        # Get comma-separated list of genes for each pathway
        path_genes = df.loc[pathway]['matching proteins in your input (labels)'].split(',')
        # Add to global set
        all_genes.update(set(path_genes))

    # Select genes with mean expression > boundary
    mask = matrix.mean(axis=0) > boundary
    strong_genes = matrix.columns[mask]

    # Keep only overlapping genes
    selected_genes = strong_genes.intersection(all_genes)
    
    return list(selected_genes)

path_genes = extract_genes(string14, pathway1, expression_matrix, 11.85)  # filtered pathway genes

# Reindex LFC matrix by numeric cluster IDs
LFC_matrix.index = [1.0,10.0,11.0,12.0,13.0,14.0,15.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
sigpath = LFC_matrix.copy()
sigpath.reset_index(inplace=True)
sigpath['index'] = sigpath['index'].astype('float64')
sigpath.sort_values(by = 'index', inplace = True)
sigpath.set_index('index', inplace = True)

# Load AZ compound cluster metadata
AZclustered = pd.read_csv('/Users/shaon/Desktop/PROTACS/github/data/AZcompound_metadata_clustered_240611a.tsv',index_col=0)

# Clean IC50 values (remove ">" then convert to float)
AZclustered['Gal'] = AZclustered['Gal'].str.replace('>','').astype(float)

# Compute cluster-wise mean IC50 and SEM
Gal = AZclustered.groupby('Dend')['Gal'].mean()
Error = AZclustered.groupby('Dend')['Gal'].sem()

# Attach IC50 means to LFC matrix
sigpath['Gal'] = Gal.values

# Plot IC50 by Cluster 
bars = pd.DataFrame({'IC50':sigpath['Gal'], 'error': Error}, index = sigpath.index)
bars.sort_values(by = 'IC50', ascending = False, inplace = True)  # sort by IC50
cmap = plt.get_cmap('jet')
bars.drop([3.0, 13.0], axis = 0, inplace = True)  # drop two clusters

# Cluster IDs for color scaling
y = list(bars.index)
rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

# Set plotting parameters
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'

# Bar plot with error bars
bars['IC50'].plot(kind = 'bar', yerr = bars['error'], capsize = 2, color = cmap(rescale(y)))
plt.title('Mitochondrial Toxicity')
plt.ylabel('IC50 (uM)') 
plt.xlabel('Cluster')
plt.savefig(outpath + 'protacs_12_barplot_IC50.pdf')
#plt.show()

# Regression modelling Proteome to IC50 
tomodel = sigpath.copy()
# tomodel.drop([13.0], inplace = True)  # optional drop

# X = pathway genes, Y = log10(IC50)
X = tomodel[path_genes]   # predictor variables
X = sm.add_constant(X)    # add intercept
Y = np.log10(tomodel['Gal'])  # response variable

gene_models = {}
p_values = {}

# Fit OLS regression for each gene individually
for gene in X.columns:  
    X_gene = X[[gene]]                # single gene as predictor
    X_gene = sm.add_constant(X_gene)  # add intercept
    model = sm.OLS(Y, X_gene)         # build model
    result = model.fit()              # fit model
    
    gene_models[gene] = result        # store fitted model
    p_values[gene] = result.pvalues[gene]  # store p-value

# Evaluate R^2 and RMSE for each gene model
r2_scores = []
rmse_scores = []

for gene, result in gene_models.items():
    X_gene = X[[gene]]
    Y = Y
    Y_pred = result.predict(sm.add_constant(X_gene))  # model predictions
    
    r2 = result.rsquared  # R² score
    rmse = np.sqrt(np.mean((Y - Y_pred) ** 2))  # RMSE
    
    r2_scores.append(r2)
    rmse_scores.append(rmse)

# Pair gene names with performance metrics
gene_r2_tuples = list(zip(gene_models.keys(), r2_scores))
gene_rmse_tuples = list(zip(gene_models.keys(), rmse_scores))

# Sort by R^2 (abs value) and RMSE
sorted_gene_r2 = sorted(gene_r2_tuples, key=lambda x: abs(x[1]))
sorted_gene_rmse = sorted(gene_rmse_tuples, key=lambda x: x[1])

# Unzip into sorted lists
sorted_gene_names_r2, sorted_r2_scores = zip(*sorted_gene_r2)
sorted_gene_names_rmse, sorted_rmse_scores = zip(*sorted_gene_rmse)

# Plot sorted R^2 scores 
plt.figure(figsize=(10, 10))
plt.plot(range(len(sorted_r2_scores)), sorted_r2_scores, marker='o', linestyle='-', color='b')
plt.title(str(pathway1[0]) +  '\n' + str(len(sorted_r2_scores)) + ' Proteins Fitted to Mitotoxicity', fontsize = 26)
plt.ylabel('$R^2$', fontsize = 26)
plt.xlabel('Sorted $R^2$ Scores', fontsize = 26)
plt.xticks(None)
plt.ylim(0,0.6)
plt.savefig(outpath + f"protacs_12_regression_models({pathway1[0].replace(' ','_')}).pdf")

sorted_gene_names_r2  # sorted list of genes by R^2

# Sort p-values by their values
sorted_genes_pvals = sorted(p_values.items(), key=lambda x: x[1])

# Extract sorted gene names and their corresponding p-values
sorted_gene_names_pval, sorted_p_values = zip(*sorted_genes_pvals)

# Convert p-values to -log10(p-value) for better visualization
log_p_values = [-np.log10(p) for p in sorted_p_values]

# Plot sorted p-values (-log10(p-values))
plt.figure(figsize=(10, 10))
plt.plot(range(len(log_p_values)), log_p_values, marker='o', linestyle='-', color='g')
plt.title('Sorted P-values for Genes', fontsize=20)
plt.ylabel('-log10(P-value)', fontsize=16)
plt.xlabel('Sorted Genes', fontsize=16)
plt.xticks(None)
plt.savefig(outpath + 'protacs_12_regression_p_values.pdf')

# Return sorted p-values for further use
sorted_gene_names_pval, sorted_p_values

# %% P-value and R^2 sorted df
df_p = pd.DataFrame({'pval': sorted_p_values}, index = sorted_gene_names_pval)
df_r = pd.DataFrame({'r2': sorted_r2_scores}, index = sorted_gene_names_r2)
df_combined = pd.merge(df_p, df_r, left_index=True, right_index=True)
df_combined.loc[['MACROH2A1','MCM6','RFC4','NSUN4','RPL11','NDUFA5']]

# %% Choose a gene for visualization of measured vs predicted pIC50 values using regression model

# Gene select (only last assignment takes effect)
#selected_gene = 'RFC4'
selected_gene = 'MCM6'
selected_gene = 'NSUN4'
selected_gene = 'NDUFA5'
#selected_gene = 'MACROH2A1'
#selected_gene = 'RPL11'

# Extract predictor variable (expression of selected gene) and response (pIC50)
X_gene = X[[selected_gene]]
Y = Y

# Add constant (intercept) to predictors
X_gene = sm.add_constant(X_gene)

# Fit linear regression model
model = sm.OLS(Y, X_gene)
result = model.fit()

# Predict pIC50 from fitted model
Y_pred = result.predict(X_gene)

# Define colors for chemical series (not used in final plots, but stored)
color = ['red', 'darkorange', 'gold', 'gold',  'gold', 'gold','gold', 'gold', 'gold', 
         'khaki', 'red', 'greenyellow', 'green', 'lightseagreen', 'blue']
color.reverse()

# Plot observed vs predicted pIC50 for the selected gene
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(Y, Y_pred, s=200, c = 'grey')  # scatter of observed vs predicted
max_value = max(Y.max(), Y_pred.max())
min_value = min(Y.min(), Y_pred.min())

# Annotate each point with its index (chemical series number)
for i, (observed, predicted) in enumerate(zip(Y, Y_pred)):
    ax.text(observed, predicted, str(i+1), fontsize=10, ha='right', va='bottom')

# Add diagonal reference line (perfect prediction)
ax.plot([min_value, max_value], [min_value, max_value], 'k--', c = 'blue')

# Title and axis labels
plt.title('Linear modelling using {} fold change\n'.format(selected_gene), fontsize = 26)
plt.xlabel('Measured Mitochondrial Toxicity (pIC50)', fontsize = 26)
plt.ylabel('Predicted Mitochondrial Toxicity (pIC50)', fontsize = 26)
plt.grid(True)

# Save observed vs predicted plot
plt.savefig(outpath + 'protacs_12_regression_' + selected_gene + '.pdf')
#plt.show()

# %% Gene expression vs IC50 scatter plot
X_val = X[[selected_gene]]  # expression values of the selected gene
Y = Y

# Configure plotting styles
plt.rcParams['axes.labelsize'] = '30'   
plt.rcParams['axes.titlesize'] = '30' 
plt.rcParams['xtick.labelsize'] = '26'  
plt.rcParams['ytick.labelsize'] = '26' 
plt.rcParams['ytick.labelsize'] = '26'

# Scatter plot: gene expression vs observed pIC50
plt.figure(figsize=(10,8))
plt.scatter(X_val, Y, c = 'grey', s = 200)

# Add regression line (predicted pIC50)
plt.plot(X_val, Y_pred, color = 'grey', linewidth = 1, ls = '--')

# Title and axis labels
plt.title('OLS using {}'.format(selected_gene) + ' fold change')
plt.xlabel('LFC')      # log fold change of gene expression
plt.ylabel('pIC50')    # negative log IC50 (toxicity)
plt.grid(True)

# Save gene expression vs IC50 plot
plt.savefig(outpath + 'protacs_12_regression_' + selected_gene + '_pred.pdf')
#plt.show()

# %% DEG p-value count extraction
# adjLFC_matrix = pd.read_csv(filepath + 'Cluster_LFCxadjPval_10uM_250305a.csv',
#                      delimiter = ';', decimal = ',', index_col=0, header = 0).T
# adjLFC_matrix.index = [1.0,10.0,11.0,12.0,13.0,14.0,15.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
# Deg_mask = np.abs(adjLFC_matrix) > -np.log10(0.05)
# df_count = Deg_mask.sum(axis = 1)
# arpro_count = df_count.drop([4,7,13]).mean()

