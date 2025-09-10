# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

data_out = os.path.dirname(__file__)
figure_out = os.path.join(data_out, '..', 'figures')
data_out = os.path.join(data_out, '..', 'data')

# %% extract input

filepath = data_out + '/'

filepath2 = figure_out + '/'

def idx_clean(df):
    df.index = df.index.str.replace('Drug_','')
    df.index = df.index.str.replace(' - DMSO','')
    return df

#%% load LFC 

t_matrix = idx_clean(pd.read_csv(filepath + 'FDA_LimmaMatrix_250304a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

adj_matrix = idx_clean(pd.read_csv(filepath + 'FDA_adjLimmaMatrix_250304a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

metadata = pd.read_csv(filepath + 'FDA_LimmaMetadata_250304a.csv', delimiter = ';', index_col = 0)

# %% Find sig counts

def counter(drug, df, cutoff):

    mask = abs(df) > cutoff
    count = mask.sum(axis = 1)
    ids = count > 1
    diffexp = ids.sum()
    diffcount_avg = count.loc[ids].mean()
    drugcount = count.loc[drug]
    drugpercent = drugcount/df.shape[1] * 100
    print('sig count, average deg, ' + drug + ' deg, percent: ')
    return diffexp, diffcount_avg, drugcount, drugpercent

counter('Methotrexate', adj_matrix, -np.log10(0.05))

counter('Doxorubicin..Adriamycin..HCl', adj_matrix, -np.log10(0.05))

# %% Counter table

# condition the p-values based on sig. cut-off
mask = abs(adj_matrix) > -np.log10(0.05)

# counting, normalizing, then conditioning based on 0 cut-off
counter = mask.sum(axis = 1).sort_values()/mask.shape[1]
#counter = mask.sum(axis = 1).sort_values()
conditioned = counter.loc[counter>0]

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
    plt.savefig(os.path.join(figure_out,saveout))
    plt.show()

    return discretized

pmf(counter, 'PMF_FDA_full.pdf')

counter.to_csv(os.path.join(data_out, 'FDA_proba_250304.tsv'))

pmf(conditioned, 'PMF_FDA_conditioned.pdf')

sns.kdeplot(counter, fill = True, clip = (0, None))
sns.kdeplot(counter, fill = True, clip = (0, counter.mean()), color = 'Green')
sns.kdeplot(counter, fill = True, clip = (counter.max(),None), color = 'Orange')
plt.axvline(counter.mean(), color = 'red', linestyle = '--')
plt.savefig(os.path.join(figure_out, 'PDF_FDA_stats.pdf'))


# %% Plot PDF 

# set target frame

df1 = counter.copy()
df2 = counter[counter == 0]

# Plot density for each column with specified colors
plt.figure(figsize=(8,6))  # Adjust figure size for better visualization

sns.kdeplot(df1.dropna(), clip=(0, None), fill=True, alpha=0.5, color='grey')
plt.axvline(0)
plt.axvline(df1['Clotrimazole'])

# Save the plot as a PDF
plt.savefig(filepath2 + 'PDF_all.pdf')
plt.show()


# %% Sort drug labels on t-Matrix

# create PCA plotting df with metadata

pca = PCA(n_components = 3)

LFC_pca = pca.fit_transform(t_matrix)

df_pca = pd.DataFrame({'PC1' : LFC_pca[:,0], 'PC2': LFC_pca[:,1], 'PC3': LFC_pca[:,2]},
                     index=t_matrix.index)

df_pca = pd.merge(df_pca, metadata, left_index=True, right_on=['Drug_'], how='left')

#  isolate drugs by class > 5 N

df_pca['n5'] = df_pca['Target_']

ncounts = df_pca['n5'].value_counts()

fivecounts = ncounts[ncounts > 4].index

fivecounts = ['Anti-infection', 'Topoisomerase', 'DHFR', '5-HT Receptor',  'Estrogen/progestogen Receptor', 'COX', 'PDE', ]

df_pca.loc[~df_pca['Target_'].isin(fivecounts), 'n5'] = 'Other'

df_pca['n5'] = pd.Categorical(df_pca['n5'], categories = fivecounts, ordered = True)

df_pca['colors'] = df_pca['n5'].cat.codes

#  set plotting labels

df_pca['labels'] = df_pca['Drug_']

tolabel = ['Amonafide', 'Doripenem.Hydrate', 'Methotrexate', 'Doxorubicin..Adriamycin..HCl', 'Pirarubicin', 'Epirubicin.HCl', 'Ibuprofen.', 'Clotrimazole', 'Ondansetron.HCl']

#tolabel = df_pca['Drug_']

df_pca.loc[~df_pca['Drug_'].isin(tolabel), 'labels'] = ''

df_pca['labels'] = df_pca['labels'].str.split('.').str[0]

df_pca['Drug_'] = df_pca['Drug_'].str.replace('/progestogen', 'Progesterone')

# %% 2D PCA

df_labels = df_pca.loc[df_pca['labels']!='']

df_empty = df_pca.loc[df_pca['n5'].isna()]

df_class = df_pca.loc[~df_pca['n5'].isna()]

explained_variance_ratios = pca.explained_variance_ratio_

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot()

scatter = ax.scatter(df_class['PC1'], df_class['PC2'], c=df_class['colors'], cmap='rainbow', s = 200, alpha = 1)

scatter2 = ax.scatter(df_empty['PC1'], df_empty['PC2'], c='grey', s = 100, alpha = 0.2)

# Annotate each point with the corresponding target label
for i, label in enumerate(df_labels['labels']):

    ax.text(df_labels['PC1'][i]- 5, df_labels['PC2'][i] + 1, label, fontsize=10)

ax.set_xlabel(f"PC1 ({explained_variance_ratios[0]*100:.2f}%)")

ax.set_ylabel(f"PC2 ({explained_variance_ratios[1]*100:.2f}%)")

handles, labels = scatter.legend_elements(prop="colors", alpha=1, size = 20)

labels = fivecounts

legend = ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', title = 'Drug Target', title_fontsize = 26, fontsize = 26)

#plt.title('Principal Component Analysis')

plt.show()

fig.savefig(filepath2 + '2D_PCA_FDA_tmatrix.pdf')

# %% 3D PCA

df_labels = df_pca.loc[df_pca['labels']!='']

df_empty = df_pca.loc[df_pca['n5'].isna()]

df_class = df_pca.loc[~df_pca['n5'].isna()]

explained_variance_ratios = pca.explained_variance_ratio_

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(111, projection='3d')

scatter = ax.scatter(df_class['PC1'], df_class['PC3'], df_class['PC2'], c=df_class['colors'], cmap='rainbow', s = 100, alpha = 1)

scatter2 = ax.scatter(df_empty['PC1'], df_empty['PC3'], df_empty['PC2'], c='grey', s = 100, alpha = 0.2)

# Annotate each point with the corresponding target label
for i, label in enumerate(df_labels['labels']):

    ax.text(df_labels['PC1'][i]- 5, df_labels['PC3'][i] + 1, df_labels['PC2'][i] - 3, label, fontsize=26)

ax.set_xlabel(f"PC1 ({explained_variance_ratios[0]*100:.2f}%)")

ax.set_zlabel(f"PC2 ({explained_variance_ratios[1]*100:.2f}%)")

ax.set_ylabel(f"PC3 ({explained_variance_ratios[2]*100:.2f}%)")

#ax.set_xticklabels([])

#ax.set_yticklabels([])

#ax.set_zticklabels([])

ax.grid(False)
ax.zaxis.pane.fill = False

handles, labels = scatter.legend_elements(prop="colors", alpha=1, size = 20)

labels = fivecounts

legend = ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', title = 'Drug Target', title_fontsize = 26, fontsize = 26)

plt.title('Principal Component Analysis')

plt.show()

fig.savefig(filepath2 + '3D_FDA_PCA_tmatrix.pdf')

# %% Export datasets

# create df with count info by drug
count_matrix = pd.DataFrame({'deg count': mask.sum(axis = 1),
'total proteins': mask.shape[1], 
'probability': mask.sum(axis=1)/mask.shape[1]
})

# merge with metadata
tableS1 = pd.merge(df_pca[['Drug_','Target_']], 
count_matrix, 
left_on = 'Drug_',
right_index = True,
how = 'inner')

#reset index
tableS1.index = tableS1['Drug_']
tableS1.drop('Drug_', axis = 1, inplace = True)

#rename 
tableS1.columns = ['Drug Class/Target', 'DEG count', 'Total Proteome', 'Probability']
tableS1.index.name = 'Approved Molecule'

# merge with LFC-signed -log10 adjusted p-value matrix
tableS1full = pd.merge(tableS1, adj_matrix, left_index = True, right_index = True)

#saveout
tableS1full.to_csv(os.path.join(data_out, 'NCB_ProteomeGuidedDiscovery_TableS1_250606a.csv'))