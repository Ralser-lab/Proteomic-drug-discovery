# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import gseapy as gp
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import leaves_list
import os
import sys
import re
import importlib
sys.path.append(os.path.join(os.path.dirname(__file__)))
import device_summarystatistics

data_out = os.path.dirname(__file__)
figure_out = os.path.join(data_out, '..', 'figures')
data_out = os.path.join(data_out, '..', 'data')

# %% Set Pathways

export_path = data_out + '/'

path = os.path.join(data_out, 'SB_FDA_metadata_250304a.tsv')

# %% Load dataset

pasef_summarized = pd.read_csv(os.path.join(data_out, 'SB_FDA_prmatrix_filtered_50_imputed_50_ltrfm_batched_summarized_250304.tsv'),
                     decimal='.', 
                     delimiter=',', 
                     index_col = 0)

metadata = pd.read_csv(path, index_col = 0, delimiter = ',')

# %% Reorder the heatmap (based on sns.clustermap Euclidian Ward) q

euclidian_ward = sns.clustermap(pasef_summarized)

row_order = leaves_list(euclidian_ward.dendrogram_row.linkage)

col_order = leaves_list(euclidian_ward.dendrogram_col.linkage)

pasef_summarized_clustered = pasef_summarized.iloc[row_order, col_order]

#  Extract sample data, batch data

metadata.loc[metadata['Drug'].isna(), 'Sample.Type'] = 'Control'

metadata = metadata[~metadata.index.duplicated(keep='first')]

metadata = metadata.loc[pasef_summarized_clustered.index]

metadata.loc[metadata['Content'].isna(), 'Sample.Type'] = 'Control'

samples = pd.get_dummies(metadata['Sample.Type'].fillna('Biological'))

batch = metadata['MS.Batch']

plate = metadata['Plate']

metadata_encoded = samples

metadata_encoded.columns = ['Drug', 'DMSO']

# Plot protein cluster heatmap (Euclidian + Ward Linkage, with orthogonal metadata)

# Calculate x and y axis tick placement 
shape = pasef_summarized_clustered.shape

xtick_gap = int(np.ceil(shape[1] / 50.0))

ytick_gap = int(np.ceil(shape[0] / 25.0))

# Start plotting

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['legend.fontsize'] = '28'  

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26'

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 10), gridspec_kw={'width_ratios': [10, 0.25, 0.5], 'wspace': 0.05})

# Plot the main protein expression heatmap

sns.heatmap(pasef_summarized_clustered, cmap='inferno', ax=ax1, cbar=False,
            xticklabels=xtick_gap, yticklabels=ytick_gap)

ax1.set_xticks(np.arange(0, shape[1], xtick_gap))

ax1.set_yticks(np.arange(0, shape[0], ytick_gap))

ax1.set_title('Heatmap (Euclidian + Ward)')

ax1.set_xlabel('Log2 Protein Abundance (' + str(shape[1]) + ')')

ax1.set_ylabel('Samples (' + str(shape[0]) + ')')

colors = plt.cm.inferno(np.linspace(0, 1, 5)) 

labels = ['2.0', '6.0', '10.0', '14.0', '18.0'] 

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]

ax1.legend(handles=patches, loc='lower left')

# Plot the one-hot encoded metadata next to the heatmap

sns.heatmap(batch.to_frame(), cmap = 'viridis', ax = ax2, cbar = False, yticklabels = False, xticklabels = False)

ax2.set_title('Batch', fontsize = 24)

# heat map extra plot for euclidian clustering by MSQC and sample

sns.heatmap(metadata_encoded, cmap='Greys', ax=ax3, cbar=False, yticklabels=False)

ax3.set_xticks(np.arange(len(metadata_encoded.columns)))

ax3.set_xticklabels(metadata_encoded.columns.values, rotation = 45, fontsize = 24)

plt.savefig(os.path.join(figure_out, 'heatmap_FDA.pdf'))

# %% Extract PCA loadings

output = pasef_summarized_clustered.copy()

output_scaled = (output - output.mean()) / output.std()

output = output_scaled.copy()

pca = PCA(n_components = 20)

pca_result = pca.fit_transform(output)  

explained_variance_ratios = pca.explained_variance_ratio_

# Plot Scree Plot

n_components = np.arange(len(explained_variance_ratios)) + 1

cumulative_variance = np.cumsum(explained_variance_ratios)

fig = plt.figure(figsize = (10,10))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

# Creating the scree plot
colors = ['skyblue' if i < 4 else 'grey' for i in n_components]

bar = plt.bar(n_components, explained_variance_ratios*100, color=colors)

plt.title('Scree Plot')

plt.xlabel('Principal Component')

plt.ylabel('Explained Variance Ratio (%)')

plt.xticks(n_components)

plt.legend(loc='upper left')

plt.savefig(os.path.join(figure_out, 'global_PCA_scree_FDA.pdf'))

# Show the plot
#plt.show()

#  %% PCA with drug targets

x = pca_result[:, 0]

y = pca_result[:, 1]

z = pca_result[:, 2]

df_pca = pd.DataFrame({'PC1':x, 'PC2':y, 'PC3':x, 'Target':metadata['Target'], 'Drug': metadata['Drug']}, index = pasef_summarized_clustered.index)

df_pca.loc[df_pca['Target'].isna(),'Target'] = 'DMSO'

df_pca.loc[df_pca['Drug'].isna(),'Drug'] = 'DMSO'

# Plot PCA with drug labels
 
plt.figure(figsize=(10, 10))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.scatter(x=df_pca['PC1'], y=df_pca['PC2'], c=df_pca['Drug'].astype('category').cat.codes, cmap='jet_r')

# Annotate each point with the corresponding target label
for i, label in enumerate(df_pca['Drug']):

    plt.annotate(label, (df_pca['PC1'][i], df_pca['PC2'][i]), fontsize=9)
    
plt.xlabel('PC1')

plt.ylabel('PC2')

plt.title('PCA Scatter Plot')

#plt.show()

# Plot PCA with target labels

plt.figure(figsize=(10, 10))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.scatter(x=df_pca['PC1'], y=df_pca['PC2'], c=df_pca['Target'].astype('category').cat.codes, cmap='jet_r')

# Annotate each point with the corresponding target label
for i, label in enumerate(df_pca['Target']):

    plt.annotate(label, (df_pca['PC1'][i], df_pca['PC2'][i]), fontsize=9)
    

plt.xlabel('PC1')

plt.ylabel('PC2')

plt.title('PCA Scatter Plot')

plt.savefig(os.path.join(figure_out, 'global_PCA_FDA.pdf'))

#plt.show()

# %% Dipsersion plots (Kernal Density Estimations, DMSO vs Controls)

pasef_summarized = pd.read_csv(os.path.join(data_out, 'SB_FDA_prmatrix_filtered_50_imputed_50_ltrfm_batched_summarized_250304.tsv'),
                     decimal='.', 
                     delimiter=',', 
                     index_col = 0)

importlib.reload(device_summarystatistics)

DMSO_frame = device_summarystatistics.calculate_cv(pasef_summarized[metadata_encoded['DMSO']==True], 'dmso')

drug_frame = device_summarystatistics.calculate_cv(pasef_summarized[~metadata_encoded['DMSO']==True], 'drug')

all_df = pd.concat([DMSO_frame, drug_frame]).reset_index(drop = True)

# %%

plt.figure(figsize=(12, 10))

custom_palette = ['blue','red'] 

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26' 

sns.jointplot(all_df, x = 'Means', y = 'Stdev', hue = 'ID', linewidth = 2, 
              kind = 'kde', space=0, height=10, ratio=4, palette = custom_palette)

sns.scatterplot(all_df, x = 'Means', y = 'Stdev', hue = 'ID', alpha = 0.2, 
                marker = 'x', palette = custom_palette)

plt.xlim(7, 15)

plt.ylim(-0.05, 0.45)

plt.legend(fontsize = '26')

plt.xlabel('Mean Protein Abundance (' + str(pasef_summarized_clustered.columns.size) + ')')

plt.ylabel('Standard Deviation (' + str(pasef_summarized_clustered.columns.size) + ')')

plt.savefig(os.path.join(figure_out, 'dispersion_kde_FDA.pdf'))
