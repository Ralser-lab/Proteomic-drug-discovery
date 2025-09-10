# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import numpy as np
from matplotlib.lines import Line2D    
from scipy.spatial import distance
import re
import sys
import importlib
sys.path.append(os.path.join(os.path.dirname(__file__)))
import device_summarystatistics

# Performs PCA on limma output on PROTAC dataset

# %% extract input
# Define file path
dir1 = os.path.dirname(__file__)
filepath2 = os.path.join(dir1, '..', 'data/')
filepath3 = os.path.join(dir1, '..', 'figures/')

# Load data
t_matrix = pd.read_csv(filepath2 + 'Cluster_LFCxPval_10uM_250305a.csv',
                       delimiter=';', decimal=',', index_col=0, header=0).T

# Perform PCA
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

#  Map Cluster labels

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
plt.rcParams['axes.titlesize'] = '30'
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
# Extracting the first occurrence of each label to ensure each label is represented once

df_pca['labels'] = df_pca['labels'].astype('float64')

df_pca.sort_values(by = ['PC1'], inplace = True)

label_colors = df_pca.drop_duplicates('clusterlabs').set_index('clusterlabs')['color']
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=12) for color in label_colors]
legend = ax.legend(handles, label_colors.index, bbox_to_anchor=(1.05, 1), loc='upper left', title='Drug Target', title_fontsize=26, fontsize=26)

plt.title('Principal Component Analysis')
plt.show()

fig.savefig(filepath3 + '3D_PCA_PROTACs_ChemicalSeries.pdf')
# %% Counter table

# clean and import data
def clean_drug_index(df):

    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)

    df.index = df.index.str.replace(' - Compound', '')

    return df

LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(filepath2, 'Drug_LFC_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

LFCxPval_matrix = clean_drug_index(pd.read_csv(os.path.join(filepath2, 'Drug_LFCxPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

LFCxadjPval_matrix = clean_drug_index(pd.read_csv(os.path.join(filepath2, 'Drug_LFCxadjPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

expression_matrix =  pd.read_csv(os.path.join(filepath2,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T

# clean and import metadata

AZmeta= pd.read_csv(os.path.join(filepath2, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)

AZmeta.index = AZmeta.index.str.replace('-','_')

# condition the p-values based on cutoff
mask = abs(LFCxadjPval_matrix) > -np.log10(0.05)

# counts after conditioning 
counts = mask.sum(axis = 1).sort_values()/mask.shape[1]

ontargetHBD = counts.loc[AZmeta['Drug_Type'] == 'Txn-PROTAC']

ARHBD = counts.loc[AZmeta['Drug_Type'] == 'AR-PROTAC']


# function to discretize and produce PMF weights    
def pmf(df, saveout):

    # outcome frequencies
    outcomes = df.value_counts().sort_index()

    # discretize outcomes
    intervals = pd.cut(outcomes.index, np.arange(-0.00001,df.max(), df.max()/100))

    discretized = outcomes.groupby(intervals, observed = False).sum()

    # calculate mass
    discretized = discretized / len(discretized)

    plt.stem(np.arange(0.00001, df.max(), df.max()/100), discretized, '-', markerfmt='o')
    plt.axvline(df.mean())
    plt.savefig(filepath3 + saveout)
    plt.show()

    return discretized

pmf(ontargetHBD, '/PMF_PROTACs_full.pdf')

pmf(ARHBD, '/PMF_PROTACs_conditioned.pdf')

# %%

# load diffexp probabilities from FDA drug screen dataset
counter = pd.read_csv(os.path.join(filepath2, 'FDA_proba_250304.tsv'), index_col = 0).iloc[:,0]

# Example standardized DataFrame (replace with your actual data)
df = pd.DataFrame({
    'FDA approved drugs': counter,
    'AR-directed HBDs': ARHBD,
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
plt.title("Density Plot with Custom Colors")

# Save the plot as a PDF
plt.savefig(filepath3 + 'PDF_FDA_PROTACs_Zcentered.pdf')
plt.show()

# %% plot PDF as boxplot
from matplotlib import pyplot as plt
df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df.items()]))

# Prepare data by dropping NaN values for each column
data_no_nan = [df[col].dropna().values for col in df.columns]

# Create the boxplot
plt.figure(figsize=(8, 6))
plt.violinplot(data_no_nan)
plt.boxplot(data_no_nan, labels = df.columns)
plt.title("Boxplot with Columns of Varying Lengths")
plt.ylabel("Values")
plt.savefig(os.path.join(filepath3, 'PDF_violinbox_FDA_PROTACs_all.pdf'))
plt.show()

# z-test with CDF
device_summarystatistics.z_test(df.iloc[:,0], df.iloc[:,1])

# %% Export datasets

# create df with count info by drug
count_matrix = pd.DataFrame({'deg count': mask.sum(axis = 1),
'total proteins': mask.shape[1], 
'probability': mask.sum(axis=1)/mask.shape[1]
})

#%%
# merge with metadata
tableS2 = pd.merge(AZmeta.loc[count_matrix.index][['Drug ID','Drug_Type', 'Binned_ligase','Binned_Target','Dend','Gal']], 
                count_matrix, 
                left_index = True,
                right_index = True,
                how = 'inner')


#rename 
tableS2.columns = ['Drug Name', 'Drug Type', 'Recruiter', 'Binder', 'Cluster', 'Gal IC50', 'DEG count', 'Total Proteome', 'Probability']
tableS2.replace(['AR-PROTAC','Non-PROTAC','Txn-PROTAC'], ['AR-HBD','not-HBD','ontarget-HBD'], inplace = True)

#sort cluster
tableS2.sort_values(by = 'Cluster', inplace = True)

# merge with LFC-signed -log10 adjusted p-value matrix
tableS2full = pd.merge(tableS2, LFCxadjPval_matrix, left_index = True, right_index = True)
tableS2.index.name = 'Unique Identifier'

#add Series reference
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14183816').values, 'Drug Name'] = 'Compound 1'
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14196658').values, 'Drug Name'] = 'Compound 2'
tableS2.loc[pd.Series(tableS2.index).str.contains('AZ14197166').values, 'Drug Name'] = 'Compound 3'

# check analogue edit
tableS2.loc[tableS2['Drug Name'].isin(['Compound 1', 'Compound 2', 'Compound 3'])]

#saveout
tableS2full.to_csv(os.path.join(filepath2, 'NCB_ProteomeGuidedDiscovery_TableS2_250606a.csv'))

