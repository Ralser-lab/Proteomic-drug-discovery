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
import matplotlib.patches as mpatches
import scanpy as sc
import seaborn as sns

class Batcher(Device): 

    def __init__ (self, data, path, batchID, logtransform):

        self.raw = data
        
        self.metapath = path

        self.batchID = batchID

        self.output = None

        self.batchdata = None

        if logtransform == True: 

            self.input = np.log2(self.raw)

        else: 

            self.input = self.raw.copy()

        return

    def batch_correct(self, toPlot = True):
    
        input = self.input.copy()

        metadata = pd.read_excel(self.metapath, index_col=0)

        batchdata = pd.DataFrame(metadata[self.batchID], index = metadata.index)

        batchdata = batchdata[~batchdata.index.duplicated(keep='first')]

        input_batchdata = pd.DataFrame(index = input.index)

        batchdata = pd.merge(input_batchdata, batchdata, how = 'inner', left_index= True, right_index = True)

        self.batchdata = batchdata

        adata = sc.AnnData(X=input.values, 
                   obs=batchdata, 
                   var=pd.DataFrame(index=range(len(input.columns))))

        sc.pp.combat(adata, self.batchID)

        output = pd.DataFrame(adata.X, index=adata.obs_names, columns=input.columns, dtype = float)  

        self.output = output

        if toPlot == True:

                self.CV_plots(self.input, title = 'before Combat')

                self.CV_plots(self.output, title = 'after Combat')
        
        return self.output      
    
    def CV_plots(self, df, title):
        
        input = df.copy()

        input['batchdata'] = self.batchdata[self.batchID]

        grouped = input.groupby('batchdata')

        cv_by_batch = grouped.apply(self.calculate_cv)

        cv_by_batch.drop('batchdata', axis = 1, inplace = True)

        shape = cv_by_batch.shape

        xtick_gap = int(self.round_up(shape[1]/50,-1))

        plt.figure(figsize=(20, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        ax = sns.heatmap(cv_by_batch, cbar = False, fmt=".2f", cmap='viridis')

        xticks = range(0, shape[1], xtick_gap)  

        ax.set_xticks(xticks)  # Set x-ticks positions

        ax.set_xticklabels(xticks, rotation=45) 

        colors = plt.cm.viridis(np.linspace(0, 1, 6))  # 6 discrete colors from the Viridis colormap

        labels = ['0%', '7%', '14%', '21%', '28%','35%']  # Customize based on your data range and preference

        patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]

        plt.ylabel('Batch')

        plt.title('Coefficient of Variation of Precursors by Batch ' + title)

        plt.legend(handles=patches, loc='lower right')

        plt.xlabel('Precursor')

        plt.show()

        long_df = cv_by_batch.reset_index().melt(id_vars='batchdata', var_name='Precursor', value_name='CV')

        plt.figure(figsize=(20,10))  # Adjust the size as needed

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        sns.boxplot(x='batchdata', y='CV', data=long_df, showfliers = False)

        plt.xticks(rotation=45)  # Rotate the precursor names for better readability

        plt.title('Distribution of CVs for Each Precursor Across Batches ' + title)

        plt.xlabel('Batch')

        plt.ylabel('Coefficient of Variation (%)')

        plt.show()

        return
 