# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import sys
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import leaves_list
import os
import re
import importlib
sys.path.append(os.path.join(os.path.dirname(__file__)))
import device_summarystatistics

# Assesses Biological Variation in protein expression matrix after pre-processing

# %% Set directory

dir_main = os.path.dirname(__file__)

export_path = os.path.join(dir_main, '..', 'data')

figure_out = os.path.join(dir_main, '..', 'figures')

filepath3 = str(export_path)

# %% Load dataset

pasef_summarized = pd.read_csv(os.path.join(filepath3, 'SB_PROTACs_prmatrix_filtered_5_imputed_50_ltrfm_batched_summarized_240314.tsv'),
                     decimal=',', 
                     delimiter=';', 
                     index_col = 0) 


# %% Reorder the heatmap (based on sns.clustermap Euclidian Ward)

euclidian_ward = sns.clustermap(pasef_summarized)

row_order = leaves_list(euclidian_ward.dendrogram_row.linkage)

col_order = leaves_list(euclidian_ward.dendrogram_col.linkage)

pasef_summarized_clustered = pasef_summarized.iloc[row_order, col_order]

metadata = pd.read_csv(os.path.join(export_path, 'SB_PROTACs_metadata_240611a.tsv'), index_col = 0)

# %% Extract drug cluster orthog data

drugcluster = metadata['Cluster'].fillna(value = 0)

# Extract Concentration data, batch data

metadata = metadata[~metadata.index.duplicated(keep='first')]

metadata = metadata.loc[pasef_summarized_clustered.index]

metadata = metadata[['Berlin_Sample.Type','Berlin_MS.Batch','AZ_Stock Concentration (mM)']]

samples = pd.get_dummies(metadata['Berlin_Sample.Type'].fillna('biological'))

concentration = pd.get_dummies(metadata['AZ_Stock Concentration (mM)'])

batch = metadata['Berlin_MS.Batch']

metadata_encoded = pd.merge(samples, concentration, left_index=True, right_index=True)

metadata_encoded.drop('biological', axis = 1, inplace = True)

metadata_encoded.columns = ['dmso\n', '0.1','1.0','10.0']

# %% Plot protein cluster heatmap (Euclidian + Ward Linkage, with orthogonal metadata)

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

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 10), gridspec_kw={'width_ratios': [20, 2, 0.5], 'wspace': 0.1})

# Plot the main protein expression heatmap

sns.heatmap(pasef_summarized_clustered, cmap='inferno', ax=ax1, cbar=False,
            xticklabels=xtick_gap, yticklabels=ytick_gap)

ax1.set_xticks(np.arange(0, shape[1], xtick_gap))

ax1.set_yticks(np.arange(0, shape[0], ytick_gap))

ax1.set_title('Protein Cluster Heatmap (Euclidian + Ward Linkage)')

ax1.set_xlabel('Log2 Protein Abundance (' + str(shape[1]) + ')')

ax1.set_ylabel('Samples (' + str(shape[0]) + ')')

colors = plt.cm.inferno(np.linspace(0, 1, 5)) 

labels = ['2.0', '6.0', '10.0', '14.0', '18.0'] 

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]

ax1.legend(handles=patches, loc='lower left')

# Plot the one-hot encoded metadata next to the heatmap

sns.heatmap(metadata_encoded, cmap='Greys', ax=ax2, cbar=False, yticklabels=False)

ax2.set_title('μmol')

ax2.set_xticks(np.arange(len(metadata_encoded.columns)))

ax2.set_xticklabels(metadata_encoded.columns.values, rotation = 45)

sns.heatmap(batch.to_frame(), cmap = 'jet_r', ax = ax3, cbar = False, yticklabels = False, xticklabels = False)

ax3.set_title('Batch')

plt.savefig(os.path.join(figure_out, 'heatmap_PROTACs.pdf'))

plt.show()

# %% Dispersion plot of summarized proteins

DMSO_frame = device_summarystatistics.calculate_cv(pasef_summarized_clustered[metadata_encoded['dmso\n']==True], 'dmso')

drug_frame = device_summarystatistics.calculate_cv(pasef_summarized_clustered[~metadata_encoded['dmso\n']==True], 'drug')

all = pd.concat([DMSO_frame, drug_frame]).reset_index(drop = True)

plt.figure(figsize=(12, 10))

custom_palette = ['blue','red'] 

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26' 

sns.jointplot(all, x = 'Means', y = 'Stdev', hue = 'ID', linewidth = 2, 
              kind = 'kde', space=0, height=10, ratio=4, palette = custom_palette)

sns.scatterplot(all, x = 'Means', y = 'Stdev', hue = 'ID', alpha = 0.2, 
                marker = 'x', palette = custom_palette)

plt.xlim(6, 18)

plt.ylim(-0.1, 1.0)

plt.legend(fontsize = '26')

plt.xlabel('Mean Protein Abundance (' + str(pasef_summarized_clustered.columns.size) + ')')

plt.ylabel('Standard Deviation (' + str(pasef_summarized_clustered.columns.size) + ')')

plt.savefig(os.path.join(figureout, 'dispersion_kde_PROTACs.pdf'))

# %% Extract PCA loadings

output = pasef_summarized_clustered.copy()

output_scaled = (output - output.mean()) / output.std()

output = output_scaled.copy()

pca = PCA(n_components = 20)

pca_result = pca.fit_transform(output)  

explained_variance_ratios = pca.explained_variance_ratio_

# %% Plot Scree Plot

n_components = np.arange(len(explained_variance_ratios)) + 1

cumulative_variance = np.cumsum(explained_variance_ratios)

fig = plt.figure(figsize = (16,12))

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

plt.savefig(os.path.join(figure_out, 'global_PCA_scree_PROTACs.pdf'))

# Show the plot
plt.show()

# %% Plot 3D PCA with drug concentration

x = pca_result[:, 0]

y = pca_result[:, 1]

z = pca_result[:, 2]

pca_metadata = pd.DataFrame({'Values' : 'PQC'}, index = pasef_summarized_clustered.index, dtype = 'string')

pca_metadata.loc[metadata_encoded['dmso\n'] == True, 'Values'] =  'dmso'

pca_metadata.loc[metadata_encoded['0.1'] == True, 'Values'] =  '0.1'

pca_metadata.loc[metadata_encoded['1.0'] == True, 'Values'] =  '1.0'

pca_metadata.loc[metadata_encoded['10.0'] == True, 'Values'] =  '10.0'

color_map = {
    'dmso': 'blue',
    '0.1': 'grey',
    '1.0': 'grey',
    '10.0': 'red'
}

default_color = 'lightgrey'  # For any values not in color_map

# Map each value in 'Values' to a color, using 'grey' as default for unmatched values
color_values = [color_map.get(val, default_color) for val in pca_metadata['Values']]

legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]

fig = plt.figure(figsize = (12,12))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26' 

ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, c = color_values, s = 50, alpha = 0.5)

ax.set_zlim(-40,60)
# Labeling the axes with the variance explained

ax.set_xlabel(f"\nPC1 ({explained_variance_ratios[0]*100:.2f}%)")

ax.set_ylabel(f"\nPC2 ({explained_variance_ratios[1]*100:.2f}%)")

ax.set_zlabel(f"\nPC3 ({explained_variance_ratios[2]*100:.2f}%)")

ax.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left')

# Title for the plot
plt.title("Principal Component Analysis")

plt.savefig(os.path.join(figure_out, 'global_PCA_concentration_PROTACs.pdf'))

plt.show()

# %% Plot 3D PCA with batch effect

pca_metadata = pd.DataFrame({'Values' : batch}, index = pasef_summarized_clustered.index, dtype = 'string')

label_mappings = {
    1: 'Batch 1',
    2: 'Batch 2',
    3: 'Batch 3',
    4: 'Batch 4',
    5: 'Batch 5',
    6: 'Batch 6',
    7: 'Batch 7',
    8: 'Batch 8',
    9: 'Batch 9',
    10: 'Batch 10',
    11: 'Batch 11'
}

# Prepare colors - assuming pca_metadata['Values'] are integers 0 through 10
cmap = plt.get_cmap('jet_r')

colors = [cmap(i / 11) for i in range(11)]  # Generate colors for each label

# Create legend patches
legend_patches = [mpatches.Patch(color=colors[i-1], label=label) for i, label in label_mappings.items()]

fig = plt.figure(figsize = (12,12))

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

ax = fig.add_subplot(111, projection='3d')

# Use the first three PCA components for x, y, and z axes
x = pca_result[:, 0]

y = pca_result[:, 1]

z = pca_result[:, 2]

ax.scatter(x, y, z, c = pca_metadata['Values'].astype('int'), cmap = 'jet_r', s = 50, 
           alpha = 0.5)

ax.set_zlim(-40,60)
# Labeling the axes with the variance explained

ax.set_xlabel(f"\nPC1 ({explained_variance_ratios[0]*100:.2f}%)")

ax.set_ylabel(f"\nPC2 ({explained_variance_ratios[1]*100:.2f}%)")

ax.set_zlabel(f"\nPC3 ({explained_variance_ratios[2]*100:.2f}%)")

ax.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left')

# Title for the plot
plt.title("Principal Component Analysis")

plt.savefig(os.path.join(figure_out, 'global_PCA_MSBatch_PROTACs.pdf'))

plt.show()

# %% Extract PC loadings

loadings = pca.components_

pc1_loadings = pd.Series(loadings[0, :], index=output.columns)

pc2_loadings = pd.Series(loadings[1, :], index=output.columns)

pc3_loadings = pd.Series(loadings[2, :], index=output.columns)


# %% Perform GSEA on PC loadings

def run_gsea(x):
  
  gsea_results = gp.prerank(rnk=x, 
                          gene_sets=os.path.join(dir_main, '..', 'data', 'c5.all.v2023.1.Hs.symbols.gmt'), 
                          outdir=export_path,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)
  
  gsea_df = gsea_results.res2d

  return gsea_df

load1 = run_gsea(pc1_loadings)

load2 = run_gsea(pc2_loadings)

load3 = run_gsea(pc3_loadings)

# %% Save GSEA results

load1.to_csv(filepath3 + '/PC1_GSEA_C5.csv')

load2.to_csv(filepath3 + '/PC2_GSEA_C5.csv')

load3.to_csv(filepath3 + '/PC3_GSEA_C5.csv')

# %% Merge PC1, PC2 and PC3 GSEA DFs

cutoff = 0.01

top_pathways = 3

def gsea_filter(x):

    filtered = x[x['FDR q-val'] < cutoff]


    forward = filtered.sort_values(by = 'NES')

    forward = forward.head(top_pathways)

    reverse = filtered.sort_values(by = 'NES', ascending = False)

    reverse = reverse.head(top_pathways)

    combined = pd.concat([forward, reverse])

    return combined

PC1_gsea = gsea_filter(pd.read_csv(filepath3 + '/PC1_GSEA_C5.csv', index_col = 'Term'))

PC2_gsea = gsea_filter(pd.read_csv(filepath3 + '/PC2_GSEA_C5.csv', index_col = 'Term'))

PC3_gsea = gsea_filter(pd.read_csv(filepath3 + '/PC3_GSEA_C5.csv', index_col = 'Term'))

gsea_merged1 = pd.merge(PC1_gsea, PC2_gsea, left_index = True, right_index = True, how = 'outer')

gsea_merged = pd.merge(gsea_merged1, PC3_gsea, left_index = True, right_index = True, how = 'outer')

gsea_redux = gsea_merged.iloc[:,[3,5,13,15,23,25]]

#gsea_redux = gsea_merged.iloc[:,[23,25,13,15,3,5]]

gsea_redux.columns = ['PC1_NES', 'PC1_FDR', 
                      'PC2_NES', 'PC2_FDR',
                      'PC3_NES', 'PC3_FDR']

top20_redux = gsea_redux.copy()

top20_redux.to_csv(filepath3 + 'PCA_GSEA_C5.csv')

top20_redux

# %% Plot NES of top 6 enriched pathways (PC1, PC2, PC3)

df = top20_redux

# Filter out non-numeric rows in FDR columns
df_numeric = df[df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].applymap(np.isreal)]

# Identify the minimum non-zero FDR value across the numeric rows
min_fdr = df_numeric[df_numeric > 0][['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].min().min()

# Replace 0 values in FDR columns with the identified minimum value
df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']] = df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].replace(0, min_fdr)

df = df.sort_values(by = ['PC1_NES', 'PC2_NES', 'PC3_NES'])

pcs = ['PC3', 'PC2', 'PC1']

# Adjusting the x-axis limits and positions to reduce white space between principal components
x_positions = np.linspace(0, 1, len(pcs))  # Equally spaced positions between 0 and 1

x_positions = [1.,0.5,0.]

# Setting up the figure

fig, ax = plt.subplots(figsize=(18, len(df.index) * 1.2))

# Looping through each principal component to plot
for i, (pc, pos) in enumerate(zip(pcs, x_positions)):

    nes_column = f"{pc}_NES"

    fdr_column = f"{pc}_FDR"
  
    # Filter out NaN values
    subset_df = df.dropna(subset=[nes_column, fdr_column])


    # edit pathway labels to be readable
    subset_df.index = [idx.split('_', 1)[1] if '_' in idx else idx for idx in subset_df.index]

    subset_df.index = ['_'.join(idx.split('_')[-6:]) if idx.count('_') >= 6 else idx for idx in subset_df.index]

    subset_df.index = subset_df.index.str.lower().str.replace('_', ' ')

    # Convert FDR to sizes for dots
    subset_df['Size'] = -np.log10(subset_df[fdr_column]) * 200
    
    # Plotting using circle outlines at the adjusted x positions
    scatter = ax.scatter([pos]*len(subset_df), subset_df.index, 
                         s=subset_df['Size'], c=subset_df[nes_column], 
                         cmap='bwr', vmin=-3, vmax=3, 
                         edgecolors='k', linewidths=0.1, facecolors='none', marker='o')

# Aesthetics and labels
ax.set_xticks(x_positions)

ax.set_xticklabels(pcs, fontsize = 30)

ax.set_xlim(-0.2, 1.2)  # Adjusted x-axis limits

ax.set_xlabel(None)

ax.set_ylabel(None)

ax.set_title(f'GSEA on PC loadings\nTop {top_pathways} Up Down by NES (FDR < {cutoff})')

plt.colorbar(scatter, ax=ax, label='Normalized Enrichment Score')

plt.tight_layout()

plt.savefig(os.path.join(figure_out, 'global_PCA_GSEA_PROTACs.pdf'))

plt.show()

# %% Calculate GSEA curves on PC loadings

gseaPC1 = gp.prerank(rnk=pc1_loadings, 
                          gene_sets=os.path.join(dir_main, '..', 'data', 'c5.all.v2023.1.Hs.symbols.gmt'), 
                          outdir=export_path,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)


gseaPC2 = gp.prerank(rnk=pc2_loadings, 
                          gene_sets=os.path.join(dir_main, '..', 'data', 'c5.all.v2023.1.Hs.symbols.gmt'), 
                          outdir=export_path,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)

gseaPC3 = gp.prerank(rnk=pc3_loadings, 
                          gene_sets=os.path.join(dir_main, '..', 'data', 'c5.all.v2023.1.Hs.symbols.gmt'), 
                          outdir=export_path,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)

# %% Plot curves of top 5 enriched pathways (PC1, PC2, PC3) 

#PC1
pc1_targets = ["GOBP_CELL_CELL_RECOGNITION",
                                "GOCC_TRANSLATION_PREINITIATION_COMPLEX",
                                "GOBP_BRANCHED_CHAIN_AMINO_ACID_METABOLIC_PROCESS",
                                "GOCC_ORGANELLAR_RIBOSOME",
                                "HP_HYPERAMMONEMIA"
                                ]



pc1_targets = ["GOBP_DNA_REPLICATION_INITIATION",
                                "GOBP_TOXIN_TRANSPORT",
                                "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING",
                                "GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX",
                                "HP_ABNORMALITY_OF_THE_MITOCHONDRION"
                                ]

pc1_targets.reverse()

axs = gseaPC1.plot(
   terms = pc1_targets, show_ranking=True, 
                                legend_kws={'loc': (0, -1.15),
                                            'fontsize': 12}
                                ) 

fig = plt.gcf()

plt.title('GSEA on PC1 loadings', size = 18)

fig.set_size_inches(5,8)

plt.tick_params(axis='y', labelsize=14) 

plt.savefig(os.path.join(figure_out, 'global_PC1_GSEA_PROTACs.pdf'))

plt.show()

#PC2
pc2_targets = ["GOBP_TRICARBOXYLIC_ACID_CYCLE",
                                "GOCC_MITOCHONDRIAL_MATRIX",
                                "GOBP_REGULATION_OF_DNA_TEMPLATED_DNA_REPLICATION",
                                "GOCC_NUCLEAR_REPLICATION_FORK",
                                "GOBP_DNA_REPLICATION_INITIATION"
                                ]

pc2_targets.reverse()

axs = gseaPC2.plot(terms = pc2_targets, show_ranking=True, 
                                legend_kws={'loc': (0, -1.15),
                                            'fontsize': 12}
                                ) 

fig = plt.gcf()

plt.title('GSEA on PC2 loadings', size = 18)

fig.set_size_inches(5,8)

plt.tick_params(axis='y', labelsize=14) 

plt.savefig(os.path.join(figure_out, 'global_PC2_GSEA_PROTACs.pdf'))

plt.show()

#PC3
pc3_targets = ["GOCC_ORGANELLE_INNER_MEMBRANE",
                                "GOCC_MITOCHONDRIAL_ENVELOPE",
                                "GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX",
                                "GOCC_SPLICEOSOMAL_COMPLEX",
                                "GOBP_RNA_SPLICING",
                                ]

terms = gseaPC3.res2d.Term

axs = gseaPC3.plot(terms = pc3_targets, show_ranking=True, 
                                legend_kws={'loc': (0, -1.15),
                                            'fontsize': 12}
                                ) 

fig = plt.gcf()

plt.title('GSEA on PC3 loadings', size = 18)

fig.set_size_inches(5,8)

plt.tick_params(axis='y', labelsize=14) 

plt.savefig(os.path.join(figure_out, 'global_PC3_GSEA_PROTACs.pdf'))

plt.show()


