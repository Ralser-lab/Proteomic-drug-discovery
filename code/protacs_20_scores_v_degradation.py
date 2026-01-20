#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_20_scores_vs_degradation.py
Description:
    Analyze model-derived toxicity scores across chemical series
    (clusters), identify bimodal series, export “toxic” and “safe” subsets,
    and generate summary visualizations. Also integrates wet-lab AR
    degradation IC50s to relate efficacy vs. predicted toxicity.

Author: Shaon Basu
Date: 2025-09-30

Inputs
------
- data/NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv
- data/Figure4_ARdegdata_SafetyScores_250523a.xlsx

Outputs
-------
- figures/protacs_20_joint_toxscores.pdf          # KDE of Toxic Probability vs Cluster
- figures/protacs_20_bimodal_toxscores.pdf        # Violin plot for bimodal series
- figures/protacs_20_barplot_toxscores.pdf        # Bar plot of tox scores by chemistry & category
- figures/rotacs_20_scatterplot_toxscores.pdf     # AR_DEG_IC50 vs Toxic Probability scatter
- data/HTMSdrugsafety_Toxic_Subset_250522a.csv      
- data/HTMSdrugsafety_Safe_SuperStrict_250522a.csv   
- data/HTMSdrugsafety_Safe_Strict_250522a.csv        
- data/HTMSdrugsafety_VHLsafe_250522a.csv            

Requirements
------------
Python >= 3.8
Packages: pandas, numpy, matplotlib, seaborn

"""
# %% Import packages
import os
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns

# Set relative paths
data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
fig_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')

# Load dataset
tableS3 = pd.read_csv(os.path.join(data_dir, 'NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv'))

# %% KDE Plot: Toxic Probability by Cluster
# Chemical series mapping dictionary
cluster_mappings = {1 : 'AR-VHL-Other',
               2: 'AR-T6N-Indole',
               3: 'AR-Other-Other',
               4: 'Txn-VHL/6N-Other',
               5: 'AR-VHL-Indole',
               6: 'AR-T5N-Other',
               7: 'Txn-Other-Other',
               8: 'AR-T6N/VHL-Pip',
               9: 'AR-T6N-Other',
               10: 'AR-T5N-Pip',
               11: 'AR-DHU-Pip',
               12: 'AR-Other-Pip',
               13: 'Non PROTAC',
               14: 'AR-L5N-Other',
               15: 'AR-L5N-Pip'} 

# Plot and save join probability density function of toxicity scores
sns.kdeplot(data=tableS3, x='Cluster', y='Toxic Probability', common_norm=False)
plt.title('KDE: Toxic Probability by Cluster')
plt.xlabel('Chemical Series')
plt.savefig(os.path.join(fig_dir, 'protacs_20_joint_toxscores.pdf'))
# plt.show()
plt.close()

# %% Violin plots for chemical series that display bimodality
# Group toxicity scores by chemical series
cluster_summary = tableS3.groupby('Cluster')['Toxic Probability'].describe().sort_values(
    ['mean', 'std'], ascending=False)
cluster_summary['chemistry'] = cluster_summary.index.map(cluster_mappings)
print(cluster_summary)

# Subset clusters with average tox score > troglitazone safety thresh (toxicity > 50)
cluster_filter1_idx = cluster_summary.loc[cluster_summary['mean' \
'']>0.50].index
cluster_filter1 = cluster_summary.loc[cluster_filter1_idx]

# Select clusters with minimum sample size thresh (n > 30)
n_thresh = 5
cluster_filter2_idx = cluster_filter1.loc[cluster_filter1['count']>n_thresh].index
cluster_filter2 = cluster_filter1.loc[cluster_filter2_idx]
toplot = tableS3.loc[tableS3['Cluster'].isin(cluster_filter2_idx)].copy()
toplot['Chemistry'] = toplot['Cluster'].map(cluster_mappings)

# Plot and save violin plots on toxicity scores (chemotypic, bimodal series)
# Drop series 6 (not bimodal) 
toplot_bimodal = toplot.loc[toplot['Cluster'].isin([1,10,11,14,15])]
sns.violinplot(toplot_bimodal, x = 'Chemistry', y = 'Toxic Probability',
               hue = 'Chemistry', palette = 'pastel', common_norm=False)
plt.legend().remove()
plt.xticks(rotation = 90)
plt.xlabel(None)
plt.savefig(os.path.join(fig_dir, 'protacs_20_bimodal_toxscores.pdf'))
# plt.show()
plt.close()

# Categorize into toxic and non-toxic groups / series using tox score cut-offs
export_columns = ['Unnamed: 0', 'Toxic Probability', 'Drug', 'Cluster','Chemistry', 'NDUFA5']
export_df = toplot[toplot['Chemistry']!='AR-Other-Pip']
toxic_subset = export_df.loc[export_df['Toxic Probability']>0.66].loc[export_df['NDUFA5']>0.44][export_columns]
toxic_subset.set_index('Unnamed: 0', inplace = True)
nontoxic_subset = export_df.loc[export_df['Toxic Probability']<0.35].loc[export_df['NDUFA5']<0.01][export_columns]
nontoxic_subset.set_index('Unnamed: 0', inplace = True)
VHLsubset = tableS3.loc[tableS3['Cluster']==1]
VHLexport = VHLsubset[['Unnamed: 0', 'Toxic Probability', 'Drug', 'Cluster','NDUFA5']]
VHLexport.set_index('Unnamed: 0', inplace = True)
VHLexport = VHLexport.assign(Chemistry='AR-VHL-Other').loc[VHLexport['NDUFA5']<0.1]

# Create summary dataframe
toxic_subset['Pred'] = 'Toxic'
nontoxic_subset['Pred'] = 'Safe'
VHLexport['Pred'] = 'Safe'
total_subset = pd.concat([toxic_subset, nontoxic_subset, VHLexport])

# Barplot of bimodal series toxicity scores (by chemical series and category)
sns.barplot(total_subset.loc[total_subset['Cluster'].isin([1,6,11,14,15])], x = 'Chemistry', y = 'Toxic Probability', hue = 'Pred')
plt.xticks(rotation = 90)
plt.ylim([-0.1,1])
plt.savefig(os.path.join(fig_dir, 'protacs_20_barplot_toxscores.pdf'))

# Read in wetlab source data of AR degradation efficacy (IC50s)
degrader_data = pd.read_excel(os.path.join(data_dir, 'Figure4_ARdegdata_SafetyScores_250523a.xlsx'),
                              names=['Compound', 'Unnamed: 0', 'AR_DEG_IC50', 'Label'])
degrader_data.iloc[:,1] = degrader_data.iloc[:,1].str.replace('-','_')

total_data = degrader_data.merge(total_subset, on='Unnamed: 0', how = 'left')
select_data = total_data.loc[total_data['Label'].isin(degrader_data['Label'].value_counts().index)]
select_data = total_data.loc[~total_data['Label'].isin(['Thalidomide Degrader', 'Lenalidomide Degrader'])]

# Scatter plot of toxicity scores vs AR degradation IC50s
plt.figure(figsize=(5, 5))  
sns.scatterplot(select_data, x = 'AR_DEG_IC50', y = 'Toxic Probability', hue = 'Label', s = 400)
plt.savefig(os.path.join(fig_dir, 'protacs_20_scatterplot_toxscores.pdf'))
plt.close()
