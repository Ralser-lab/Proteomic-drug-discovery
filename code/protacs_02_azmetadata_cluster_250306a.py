# %% Load packages
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster
import os

# Cleans AZ metadata (on drug level), clusters based on chemical features.
# Updates chemical cluster labels to metadata file and for dendogram plotting in R

# %% Set relative wd

wd = os.path.dirname(__file__)

data = os.path.join(wd, '..', 'data')

AZ_meta = pd.read_csv(os.path.join(data, 'AZcompound_metadata_240611a.csv'),
                header = 0, index_col = 0)

AZ_meta.drop(columns='to_explode.1', inplace=True)

# %% Create drug labels based on Drug ID and E3 Ligase Chemistry 

AZ_cluster = AZ_meta[['Drug ID', 'E3_ligase','target']].copy()

# Label non-AR-PROTAC PROTACs as 'Txn-PROTAC'

AZ_cluster.loc[AZ_cluster['Drug ID'].str.contains('PROTAC') & ~AZ_cluster['Drug ID'].str.contains('AR'), 'Drug ID'] = 'Txn-PROTAC'

# Label non-PROTACs as 'Non-PROTAC'

AZ_cluster.loc[~AZ_cluster['Drug ID'].str.contains('PROTAC'), 'Drug ID'] = 'Non-PROTAC'

# %%

# Label E3 ligase with value counts < 8 to 'Other'

value_counts = AZ_cluster['E3_ligase'].value_counts()

replace = value_counts[value_counts < 8].index

AZ_cluster.loc[AZ_cluster['E3_ligase'].isin(replace), 'E3_ligase'] = 'Other'

# Label Ligand target with value counts < 20 to 'Other'

target_counts = AZ_cluster['target'].value_counts()

replace2 = target_counts[target_counts < 20].index 

AZ_cluster.loc[AZ_cluster['target'].isin(replace2), 'target'] = 'Other'

# Label non-PROTACs to have 'None' in E3 ligase column and target columns

AZ_cluster.loc[AZ_cluster['Drug ID'] == 'Non-PROTAC', 'E3_ligase'] = 'None'

AZ_cluster.loc[AZ_cluster['Drug ID'] == 'Non-PROTAC', 'target'] = 'None'

AZ_cluster.loc[AZ_cluster['target'].isna(), 'target'] = 'Other'

# Label dihydrouracil containing PROTACs to CRBN_dihydrouracil

AZ_cluster.loc[AZ_cluster['E3_ligase'].str.contains('dihydrouracil'), 'E3_ligase'] = 'CRBN_dihydrouracil'

# Label VHL amide tBu containing PROTACs to VHL_amide_tBu

AZ_cluster.loc[AZ_cluster['E3_ligase'].str.contains('VHL'), 'E3_ligase'] = 'VHL_amide_tBu'

AZ_cluster.loc[AZ_cluster['target'].str.contains('piperidine'), 'target'] = 'AR_piperidine'

leaf_labels = AZ_cluster['Drug ID']

# %% One hot encoded drug labels

AZ_encoded = pd.get_dummies(AZ_cluster, prefix = '', prefix_sep= '')

Z = linkage(AZ_encoded, method = 'ward', metric = 'euclidean')

plt.figure(figsize=(10, 7))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Small Molecules (' + str(len(leaf_labels)) + ')')
plt.ylabel('Distance')
dendrogram(
    Z,
    labels = leaf_labels,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=2.,  # font size for the x axis labels
)
plt.show()

# %% Cut ward euclidian tree into 15 clusters

cluster_labels = fcluster(Z, 15, criterion='maxclust')

AZ_cluster['Cluster ID'] = cluster_labels

AZ_cluster['Cluster ID'].value_counts().plot(kind = 'bar')

AZ_cluster.columns = ['Drug_Type', 'Binned_ligase', 'Binned_Target', 'Cluster']

# %% Remap clusters onto PC1 projection

remapping = {
    15:1, 8:2, 14:3, 10:4, 
    7:5, 13:6, 9:7, 6:8,
    12 : 9, 5:10, 4:11, 3:12, 1:13, 11:14, 2:15
}

AZ_cluster['Dend'] = AZ_cluster['Cluster'].map(remapping)

# %% export cleaned AstraZeneca metadata 

AZ_cleaned = pd.merge(AZ_meta, AZ_cluster, left_index = True, right_index = True)

AZ_export = AZ_cleaned.copy()

AZ_export.drop(['Cluster'], axis =1, inplace = True)

AZ_export.sort_values('Dend', inplace = True)

AZ_export.rename(columns={'Drug_Type':'Drug Class','Dend':'Cluster','Binned_ligase':'Degrader','Binned_Target':'Warhead'},inplace=True)

AZ_export.index.names = ['ID']

AZ_export.loc[AZ_export['Drug Class'] == 'Txn-PROTAC', 'Drug Class'] = 'Nuclear-PROTAC'

AZ_cleaned.to_csv(os.path.join(data,'AZcompound_metadata_clustered_240611a.tsv'), index = True, header = True)

# %% export cleaned AstraZeneca metadata for circular plotting in R

AZ_forR = AZ_encoded.copy()

AZ_forR[['Drug_Type','Ligase', 'Target', 'Cluster']] = AZ_cleaned[['Drug_Type','Binned_ligase','Binned_Target','Dend']]

AZ_forR.to_csv(os.path.join(data,'AZcompound_metadata_onehotencoded_240611a.tsv'), index = True, header = True)

# %% Cluster Metaplots

figures = os.path.join(wd,'..','figures')

plt.figure(figsize = (2,2))
AZ_cluster['Drug_Type'].value_counts().plot(kind = 'bar')
plt.ylabel('Count')
plt.xlabel('')
plt.savefig(os.path.join(figures,'class_metaplot.pdf'))

AZ_metaplot = AZ_cluster[AZ_cluster['Drug_Type']!='Non-PROTAC']

AZ_metaplot = AZ_metaplot[AZ_metaplot['Binned_ligase']!='Other']

plt.figure(figsize = (2,2))
AZ_metaplot['Binned_ligase'].value_counts().plot(kind = 'bar')
plt.axhline(y=20, c = 'red', linestyle = '--', linewidth = 0.5)
plt.ylabel('Count')
plt.xlabel('')
plt.savefig(os.path.join(figures,'degrader_metaplot.pdf'))

AZ_metaplot = AZ_cluster[AZ_cluster['Binned_Target']!='Non-PROTAC']

AZ_metaplot = AZ_metaplot[AZ_metaplot['Binned_Target']!='Other']

AZ_metaplot = AZ_metaplot[AZ_metaplot['Binned_Target']!='None']

plt.figure(figsize = (2,2))
AZ_metaplot['Binned_Target'].value_counts().plot(kind = 'bar')
plt.axhline(y=20, c = 'red', linestyle = '--', linewidth = 0.5)
plt.ylabel('Count')
plt.xlabel('')
plt.savefig(os.path.join(figures,'ligand_metaplot.pdf'))
