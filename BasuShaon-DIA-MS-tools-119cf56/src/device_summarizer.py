# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device import Device

# third party imports
import gseapy as gp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import leaves_list
import seaborn as sns
from sklearn.decomposition import PCA
import subprocess

class Summarizer(Device):

    def __init__(self, df, lfqscript, directory):

        self.input = df

        self.input_long = None

        self.lfqscript = lfqscript

        self.directory = directory

        self.clustered_matrix = None

        self.format_df_maxlfq()

    def format_df_maxlfq(self):
    
        pasef_batch_long = self.input.reset_index().melt(id_vars=['index'])

        pasef_batch_long = pasef_batch_long.rename(
        columns={'index': 'sample_list', 'Precursor.Id': 'id', 'Genes': 'protein_list', 'value': 'quant'})

        pasef_batch_long.index = pasef_batch_long['id']

        pasef_batch_long.drop(axis = 1, columns = 'id', inplace = True)

        self.input_long = pasef_batch_long.copy()

        return

    def maxlfq(self, longform, convert = False):

        r_script_path = self.lfqscript

        wd = os.path.join(self.directory, 'output')

        print(wd)

        file_path = longform

        out_name = file_path.replace('long','summarized')

        # Ensure convert is passed as a string "TRUE" or "FALSE"
        convert_str = "TRUE" if convert else "FALSE"

        command = ['Rscript', r_script_path, wd, file_path, out_name, convert_str]

        process = subprocess.run(command, capture_output=True, text=True)

        print('STDOUT:', process.stdout)

        print('STDERR:', process.stderr)

        return

    def protein_cluster_matrix(self, Routput):

        protein_matrix = Routput

        shape = protein_matrix.shape

        euclidian_ward = sns.clustermap(protein_matrix)

        row_order = leaves_list(euclidian_ward.dendrogram_row.linkage)

        col_order = leaves_list(euclidian_ward.dendrogram_col.linkage)

        protein_matrix_clustered = protein_matrix.iloc[row_order, col_order]

        self.clustered_matrix = protein_matrix_clustered.copy()

        xtick_gap = int(self.round_up(shape[1]/50,-1))

        ytick_gap = int(self.round_up(shape[0]/25,-1))

        plt.figure(figsize=(20,10))  # Adjust the size as needed

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26'

        ax = sns.heatmap(protein_matrix_clustered, cmap='inferno', 
                 xticklabels=1, yticklabels=1, cbar = False)

        xticks = range(0, shape[1], xtick_gap)  

        yticks = range(0, shape[0], ytick_gap)

        ax.set_xticks(xticks)  # Set x-ticks positions

        ax.set_xticklabels(xticks, rotation=45) 

        ax.set_yticks(yticks)

        ax.set_yticklabels(yticks)

        plt.ylabel('Samples (Euclidian Ward)')

        plt.title('Protein Expression Matrix')

        plt.xlabel('Log2 Abundance (Euclidian Ward)')

        colors = plt.cm.inferno(np.linspace(0, 1, 5))  # 5 discrete colors from the Viridis colormap

        labels = ['lowest', 'low', 'medium', 'high', 'highest']  # Customize based on your data range and preference

        patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]

        plt.legend(handles=patches, loc='lower right')

        plt.show()

        self.global_variance_analysis()

        return
    
    def global_variance_analysis(self):

        output = self.clustered_matrix.copy()

        output_scaled = (output - output.mean()) / output.std()

        output = output_scaled.copy()

        pca = PCA(n_components = 20)

        pca_result = pca.fit_transform(output)  

        explained_variance_ratios = pca.explained_variance_ratio_

        n_components = np.arange(len(explained_variance_ratios)) + 1

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

        # Show the plot
        plt.show()

        x = pca_result[:, 0]

        y = pca_result[:, 1]

        z = pca_result[:, 2]
    
        fig = plt.figure(figsize = (12,12))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        plt.rcParams['ytick.labelsize'] = '26' 

        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(x, y, z, s = 50, alpha = 0.5)

        ax.set_zlim(-40,60)
        # Labeling the axes with the variance explained

        ax.set_xlabel(f"\nPC1 ({explained_variance_ratios[0]*100:.2f}%)")

        ax.set_ylabel(f"\nPC2 ({explained_variance_ratios[1]*100:.2f}%)")

        ax.set_zlabel(f"\nPC3 ({explained_variance_ratios[2]*100:.2f}%)")

        plt.title("Principal Component Analysis")

        plt.show()

        loadings = pca.components_

        self.gsea_plot(loadings)

        return

    def run_gsea(self, x, sets = '/Users/shaon/Desktop/PROTACS/PROTACS/c5.all.v2023.1.Hs.symbols.gmt', 
                 o = '/Users/shaon/Desktop/PROTACS/PROTACS/'):
  
        gsea_results = gp.prerank(rnk=x, 
                          gene_sets= sets,
                          outdir= o,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)
   
        gsea_df = gsea_results.res2d

        return gsea_df

    def gsea_filter(self, x, cutoff = 0.01, top_pathways = 3):

        filtered = x[x['FDR q-val'] < cutoff]

        forward = filtered.sort_values(by = 'NES')

        forward = forward.head(top_pathways)

        reverse = filtered.sort_values(by = 'NES', ascending = False)

        reverse = reverse.head(top_pathways)

        combined = pd.concat([forward, reverse])

        return combined

    def gsea_plot(self, loadings):

        cutoff = 0.01

        top_pathways = 3

        pcs = ['PC1', 'PC2', 'PC3']
        
        pc1_loadings = pd.Series(loadings[0, :], index=self.clustered_matrix.columns)

        pc2_loadings = pd.Series(loadings[1, :], index=self.clustered_matrix.columns)

        pc3_loadings = pd.Series(loadings[2, :], index=self.clustered_matrix.columns)

        load1_gsea = self.run_gsea(pc1_loadings)

        load2_gsea = self.run_gsea(pc2_loadings)

        load3_gsea = self.run_gsea(pc3_loadings)
        
        load1_gsea.set_index('Term', inplace = True)

        load2_gsea.set_index('Term', inplace = True)

        load3_gsea.set_index('Term', inplace = True)

        PC1_gsea = self.gsea_filter(load1_gsea, cutoff, top_pathways)

        PC2_gsea = self.gsea_filter(load2_gsea, cutoff, top_pathways)

        PC3_gsea = self.gsea_filter(load3_gsea, cutoff, top_pathways)

        gsea_merged1 = pd.merge(PC1_gsea, PC2_gsea, left_index = True, right_index = True, how = 'outer')

        gsea_merged = pd.merge(gsea_merged1, PC3_gsea, left_index = True, right_index = True, how = 'outer')

        gsea_redux = gsea_merged.iloc[:,[3,5,13,15,23,25]]

        gsea_redux.columns = ['PC1_NES', 'PC1_FDR', 
                                'PC2_NES', 'PC2_FDR',
                                'PC3_NES', 'PC3_FDR']

        # Plot NES of top 6 enriched pathways (PC1, PC2, PC3)    

        df = gsea_redux.copy()

        return df
    


        # Filter out non-numeric rows in FDR columns
        df_numeric = df[df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].applymap(np.isreal)]

        # Identify the minimum non-zero FDR value across the numeric rows
        min_fdr = df_numeric[df_numeric > 0][['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].min().min()

        # Replace 0 values in FDR columns with the identified minimum value
        df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']] = df[['PC1_FDR', 'PC2_FDR', 'PC3_FDR']].replace(0, min_fdr)

        df = df.sort_values(by = ['PC1_NES', 'PC2_NES', 'PC3_NES'])

        # Adjusting the x-axis limits and positions to reduce white space between principal components
        x_positions = np.linspace(0, 1, len(pcs))  # Equally spaced positions between 0 and 1

        # Setting up the figure

        fig, ax = plt.subplots(figsize=(18, len(df.index) * 1.2))

        # Looping through each principal component to plot
        for i, (pc, pos) in enumerate(zip(pcs, x_positions)):

            nes_column = f"{pc}_NES"

            fdr_column = f"{pc}_FDR"
  
            # Filter out NaN values
            subset_df = df.dropna(subset=[nes_column, fdr_column])

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

        plt.show()      

        return  
