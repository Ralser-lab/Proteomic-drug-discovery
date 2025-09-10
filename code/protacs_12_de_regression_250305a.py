# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# %% set directory

dir_main = os.path.dirname(__file__)

filepath = os.path.join(dir_main, '..', 'data/')

outpath = os.path.join(dir_main, '..', 'figures/')

expression_matrix =  pd.read_csv(filepath + 'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv',
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T

LFC_matrix = pd.read_csv(filepath + 'Cluster_LFC_10uM_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T
                     
string14 = pd.read_csv(filepath + 'Cluster14_enrichment.Component.tsv', index_col = 1, sep = '\t')

string4 = pd.read_csv(filepath + 'Cluster4_enrichment.Component.tsv', index_col = 1, sep = '\t')

# %%

pathway1 = list(string14.index[0:5]) + list(string4.index[0:5])

#pathway1 = ['Inner mitochondrial membrane protein complex']

#pathway1 = ['Oxidoreductase complex']

#pathway1 = ['Mitochondrial ribosome']

#pathway1 = ['Mitochondrial respirasome']

#pathway1 = ['Respiratory chain complex']

#pathway1 = ['Chromosomal region']

#pathway1 = ['Nuclear chromosome']

#pathway1 = ['MCM complex']

#pathway1 = ['Cytosolic ribosome']

#pathway1 = ['Replication fork']

def extract_genes(df, pathways, matrix, boundary):

    all_genes = set()

    for pathway in pathways:

        path_genes = df.loc[pathway]['matching proteins in your input (labels)'].split(',')

        all_genes.update(set(path_genes))

    mask = matrix.mean(axis=0) > boundary

    strong_genes = matrix.columns[mask]

    selected_genes = strong_genes.intersection(all_genes)
    
    return list(selected_genes)

path_genes = extract_genes(string14, pathway1, expression_matrix, 11.85)

# 

LFC_matrix.index = [1.0,10.0,11.0,12.0,13.0,14.0,15.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]

sigpath = LFC_matrix.copy()

sigpath.reset_index(inplace=True)

sigpath['index'] = sigpath['index'].astype('float64')

sigpath.sort_values(by = 'index', inplace = True)

sigpath.set_index('index', inplace = True)

AZclustered = pd.read_csv('/Users/shaon/Desktop/PROTACS/github/data/AZcompound_metadata_clustered_240611a.tsv',index_col=0)

AZclustered['Gal'] = AZclustered['Gal'].str.replace('>','').astype(float)

Gal = AZclustered.groupby('Dend')['Gal'].mean()

Error = AZclustered.groupby('Dend')['Gal'].sem()

sigpath['Gal'] = Gal.values

#  Plot IC50 By Cluster

bars = pd.DataFrame({'IC50':sigpath['Gal'], 'error': Error}, index = sigpath.index)

bars.sort_values(by = 'IC50', ascending = False, inplace = True)

cmap = plt.get_cmap('jet')

bars.drop([3.0, 13.0], axis = 0, inplace = True)

y = list(bars.index)

rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26'

bars['IC50'].plot(kind = 'bar', yerr = bars['error'], capsize = 2, color = cmap(rescale(y)))

plt.title('Mitochondrial Toxicity')

plt.ylabel('IC50 (uM)') 

plt.xlabel('Cluster')

plt.savefig(outpath + 'barplot_IC50.pdf')

plt.show()

# Regression modelling Proteome to IC50

tomodel = sigpath.copy()

#tomodel.drop([13.0], inplace = True)

X = tomodel[path_genes]
  # Predictor variables

X = sm.add_constant(X)

Y = np.log10(tomodel['Gal'])  # Response variable

gene_models = {}
p_values = {}

# Iterate over each gene
for gene in X.columns:  # Assuming the columns represent genes
    # Extract the predictor variables (X_gene) and response variable (Y)
    X_gene = X[[gene]]
    Y = Y
    
    # Add a constant term to the predictor variables (intercept)
    X_gene = sm.add_constant(X_gene)
    
    # Fit the linear regression model for the current gene
    model = sm.OLS(Y, X_gene)
    result = model.fit()
    
    # Store the regression results in the dictionary
    gene_models[gene] = result

    # Extract p-values and store them in the dictionary
    p_values[gene] = result.pvalues[gene]

# Create empty lists to store R^2 and RMSE for each gene
r2_scores = []
rmse_scores = []

# Iterate over each gene in gene_models
for gene, result in gene_models.items():
    # Extract predictor variables (X_gene) and response variable (Y)
    X_gene = X[[gene]]
    Y = Y
    
    # Predict the response variable using the fitted model
    Y_pred = result.predict(sm.add_constant(X_gene))
    
    # Calculate R^2 score and RMSE using the fitted model
    r2 = result.rsquared
    rmse = np.sqrt(np.mean((Y - Y_pred) ** 2))
    
    # Append scores to the lists
    r2_scores.append(r2)
    rmse_scores.append(rmse)

# Zip gene names with their corresponding R^2 scores or RMSE values
gene_r2_tuples = list(zip(gene_models.keys(), r2_scores))
gene_rmse_tuples = list(zip(gene_models.keys(), rmse_scores))

# Sort the tuples based on the second element (R^2 score or RMSE value)
sorted_gene_r2 = sorted(gene_r2_tuples, key=lambda x: abs(x[1]))
sorted_gene_rmse = sorted(gene_rmse_tuples, key=lambda x: x[1])
# Extract sorted gene names and scores
sorted_gene_names_r2, sorted_r2_scores = zip(*sorted_gene_r2)
sorted_gene_names_rmse, sorted_rmse_scores = zip(*sorted_gene_rmse)

# Plot sorted R^2 scores
plt.figure(figsize=(10, 10))
plt.plot(range(len(sorted_r2_scores)), sorted_r2_scores, marker='o', linestyle='-', color='b')
#plt.title('Network Nodes ' + str(len(sorted_r2_scores)) +  ' Fitted to \nOrthogonal Mitochondrial Toxicity', fontsize = 26)
plt.title(str(pathway1[0]) +  '\n' + str(len(sorted_r2_scores)) + ' Proteins Fitted to Mitotoxicity', fontsize = 26)
plt.ylabel('$R^2$', fontsize = 26)
plt.xlabel('Sorted $R^2$ Scores', fontsize = 26)
plt.xticks(None)
plt.ylim(0,0.6)
plt.savefig(outpath + str(pathway1[0]) + ' models.pdf')
sorted_gene_names_r2

# %%
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
plt.savefig(outpath + 'regression_p_values.pdf')

# Return sorted p-values for further use
sorted_gene_names_pval, sorted_p_values


# %% P-value and R2 sorted df

df_p = pd.DataFrame({'pval': sorted_p_values}, index = sorted_gene_names_pval)
df_r = pd.DataFrame({'r2': sorted_r2_scores}, index = sorted_gene_names_r2)
df_combined = pd.merge(df_p, df_r, left_index=True, right_index=True)

df_combined.loc[['MACROH2A1','MCM6','RFC4','NSUN4','RPL11','NDUFA5']]

# %% Choose a gene for visualization
#selected_gene = 'RFC4'
selected_gene = 'MCM6'
selected_gene = 'NSUN4'
selected_gene = 'NDUFA5'
#selected_gene = 'MACROH2A1'
#selected_gene = 'RPL11'


# Extract predictor variables (X_gene) and response variable (Y) for the selected gene
#X.index = ['Txn-PROTAC','Controls','AR-PRO VHL amide', 'AR-PRO Thalidomide 6N', 'AR-PRO Other', 'AR-PRO Lenalinomide 5N', 'AR-PRO Dihydrouracyl', 'AR-PRO Thalidomide 5N']
X_gene = X[[selected_gene]]
Y = Y
#Y.index = ['Txn-PROTAC','Controls','AR-PRO VHL amide', 'AR-PRO Thalidomide 6N', 'AR-PRO Other', 'AR-PRO Lenalinomide 5N', 'AR-PRO Dihydrouracyl', 'AR-PRO Thalidomide 5N']
# Add a constant term to the predictor variables (intercept)
X_gene = sm.add_constant(X_gene)

# Fit the linear regression model for the selected gene
model = sm.OLS(Y, X_gene)
result = model.fit()

# Predict the response variable using the fitted model
Y_pred = result.predict(X_gene)

color = ['red', 'darkorange', 'gold', 'gold',  'gold', 'gold','gold', 'gold', 'gold', 'khaki', 'red', 'greenyellow', 'green', 'lightseagreen', 'blue']  # Cluster colors for each point
color.reverse()

#Y = Y.loc[Y.index!='AR-PRO Other']

# Plot observed vs predicted values
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(Y, Y_pred, s=200, c = 'grey')
max_value = max(Y.max(), Y_pred.max())
min_value = min(Y.min(), Y_pred.min())
for i, (observed, predicted) in enumerate(zip(Y, Y_pred)):
    ax.text(observed, predicted, str(i+1), fontsize=10, ha='right', va='bottom')
ax.plot([min_value, max_value], [min_value, max_value], 'k--', c = 'blue')
plt.title('Linear modelling using {} fold change\n'.format(selected_gene), fontsize = 26)
plt.xlabel('Measured Mitochondrial Toxicity (pIC50)', fontsize = 26)
plt.ylabel('Predicted Mitochondrial Toxicity (pIC50)', fontsize = 26)
plt.grid(True)
plt.savefig(outpath + 'OLS_' + selected_gene + '.pdf')
plt.show()

rmse = np.sqrt(np.mean((Y - Y_pred) ** 2))
nrmse = rmse / (np.max(Y) - np.min(Y))

#%%
X_val = X[[selected_gene]]

Y = Y

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26'

# Plot gene expression vs IC50 values
plt.figure(figsize=(10,8))
plt.scatter(X_val, Y, c = 'grey', s = 200)
plt.plot(X_val, Y_pred, color = 'grey', linewidth = 1, ls = '--')
plt.title('OLS using {}'.format(selected_gene) + ' fold change')
plt.xlabel('LFC')
plt.ylabel('pIC50')
plt.grid(True)
plt.savefig(outpath + 'OLS_' + selected_gene + '_pred.pdf')
plt.show()

# %% DEG p-value count extraction

adjLFC_matrix = pd.read_csv(filepath + 'Cluster_LFCxadjPval_10uM_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T

adjLFC_matrix.index = [1.0,10.0,11.0,12.0,13.0,14.0,15.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]

Deg_mask = np.abs(adjLFC_matrix) > -np.log10(0.05)

df_count = Deg_mask.sum(axis = 1)

arpro_count = df_count.drop([4,7,13]).mean()

