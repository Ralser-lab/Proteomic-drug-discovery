# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import numpy as np
import gseapy as gp
import os

# Performs GSEA on limma output on PROTAC dataset
# %% extract input

dir_main = os.path.dirname(__file__)

filepath2 = os.path.join(dir_main, '..', 'data/')

filepath3 = os.path.join(dir_main, '..', 'data')

fig_out = os.path.join(dir_main, '..', 'figures')

# %% Extract 10 micromolar values
DE_1 = pd.read_csv(filepath2 + 'Cluster1_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)

DE_2 = pd.read_csv(filepath2 + 'Cluster2_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_3 = pd.read_csv(filepath2 + 'Cluster3_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_4 = pd.read_csv(filepath2 + 'Cluster4_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_5 = pd.read_csv(filepath2 + 'Cluster5_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_6 = pd.read_csv(filepath2 + 'Cluster6_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_7 = pd.read_csv(filepath2 + 'Cluster7_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)

DE_8 = pd.read_csv(filepath2 + 'Cluster8_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_9 = pd.read_csv(filepath2 + 'Cluster9_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_10 = pd.read_csv(filepath2 + 'Cluster10_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_11 = pd.read_csv(filepath2 + 'Cluster11_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_12 = pd.read_csv(filepath2 + 'Cluster12_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_13 = pd.read_csv(filepath2 + 'Cluster13_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_14 = pd.read_csv(filepath2 + 'Cluster14_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_15 = pd.read_csv(filepath2 + 'Cluster15_10uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

#  Extract 1 micromolar values 

DE_1_1 = pd.read_csv(filepath2 + 'Cluster1_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)

DE_1_2 = pd.read_csv(filepath2 + 'Cluster2_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_3 = pd.read_csv(filepath2 + 'Cluster3_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_4 = pd.read_csv(filepath2 + 'Cluster4_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_5 = pd.read_csv(filepath2 + 'Cluster5_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_6 = pd.read_csv(filepath2 + 'Cluster6_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_7 = pd.read_csv(filepath2 + 'Cluster7_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)

DE_1_8 = pd.read_csv(filepath2 + 'Cluster8_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_9 = pd.read_csv(filepath2 + 'Cluster9_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_10 = pd.read_csv(filepath2 + 'Cluster10_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_11 = pd.read_csv(filepath2 + 'Cluster11_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_12 = pd.read_csv(filepath2 + 'Cluster12_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_13 = pd.read_csv(filepath2 + 'Cluster13_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_14 = pd.read_csv(filepath2 + 'Cluster14_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_1_15 = pd.read_csv(filepath2 + 'Cluster15_1p0uM_Limma_250305a.csv',
                     delimiter = ';', decimal = ',', index_col=0)


# %% Limma DE 10 micromolar 
def run_gsea(x):
  gsea_results = gp.prerank(rnk=x, 
                          gene_sets=os.path.join(filepath3, 'c5.go.cc.v2023.2.Hs.symbols.gmt'), 
                          outdir=filepath3,
                          min_size=15,  # Minimum size of gene set to consider
                          max_size=500, # Maximum size of gene set to consider
                          processes=4)
  gsea_df = gsea_results.res2d
  return gsea_df

def gsealimma(x):
    x['sign'] = np.sign(x['logFC']) * -np.log10(x['adj.P.Val'])
    ranked = x['sign'].sort_values(ascending=False)
    #ranked = x['logFC'].sort_values(ascending=False)
    gsea_DE = run_gsea(ranked)
    return gsea_DE

#  GSEA on limma differential expression sets
df1 = gsealimma(DE_1)
df2 = gsealimma(DE_2)
df3 = gsealimma(DE_3)
df4 = gsealimma(DE_4)
df5 = gsealimma(DE_5)
df6 = gsealimma(DE_6)
df7 = gsealimma(DE_7)
df8 = gsealimma(DE_8)
df9 = gsealimma(DE_9)
df10 = gsealimma(DE_10)
df11 = gsealimma(DE_11)
df12 = gsealimma(DE_12)
df13 = gsealimma(DE_13)
df14 = gsealimma(DE_14)
df15 = gsealimma(DE_15)

# Limma DE 1 micromolar 
df16 = gsealimma(DE_1_1)
df17 = gsealimma(DE_1_2)
df18 = gsealimma(DE_1_3)
df19 = gsealimma(DE_1_4)
df20 = gsealimma(DE_1_5)
df21 = gsealimma(DE_1_6)
df22 = gsealimma(DE_1_7)
df23 = gsealimma(DE_1_8)
df24 = gsealimma(DE_1_9)
df25 = gsealimma(DE_1_10)
df26 = gsealimma(DE_1_11)
df27 = gsealimma(DE_1_12)
df28 = gsealimma(DE_1_13)
df29 = gsealimma(DE_1_14)
df30 = gsealimma(DE_1_15)

# %% Plot GSEA by NES and FDR cutoffs using 10 micromolar dataset

gsea_dict = {
    "df16": df16,
    "df17": df17,
    "df18": df18,
    "df19": df19,
    "df20": df20,
    "df21": df21,
    "df22": df22,
    "df23": df23,
    "df24": df24,
    "df25": df25,
    "df26": df26,
    "df27": df27,
    "df28": df28,
    "df29": df29,
    "df30": df30,
    "df1": df1,
    "df2": df2,
    "df3": df3,
    "df4": df4,
    "df5": df5,
    "df6": df6,
    "df7": df7,
    "df8": df8,
    "df9": df9,
    "df10": df10,
    "df11": df11,
    "df12": df12,
    "df13": df13,
    "df14": df14,
    "df15": df15}

#dogma: df containing all the gsea pathways
dogma = pd.DataFrame()

cat_order = list(gsea_dict.keys())

#cat_order.reverse()

filter = 0.001 #FDR cutoff
filter2 = 0 #NES cutoff
top = 5

#Create top GSEA table
for idx, df in gsea_dict.items():
    df_subset = df.loc[df['FDR q-val'] < filter]
    df_subset = df_subset.loc[abs(df_subset['NES']) > filter2]
    df_subset.sort_values(by='FDR q-val', ascending=True)
    df_subset = df_subset.head(top)
    df_subset.index = df_subset['Term']
    df_subset['condition'] = idx
    dogma = pd.concat([dogma,df_subset], axis =0)

dogma2 = dogma.copy()

#Split and take substring
dogma2['Term'] = [idx.split('_', 1)[1] if '_' in idx else idx for idx in dogma2['Term']]
#term filter
dogma2['Term'] = ['_'.join(idx.split('_')[-9:]) if idx.count('_') >= 6 else idx for idx in dogma2['Term']]
#lower case 
dogma2['Term'] = dogma2['Term'].str.lower().str.replace('_', ' ')
#convert NES to numeric
dogma2['NES'] = pd.to_numeric(dogma2['NES'], errors = 'coerce')
#sort df by NES
dogma2.sort_values(by = 'NES', inplace = True, ascending = False)

dogma2.to_csv(os.path.join(filepath3, 'gsea_FDR_NES_redux_HBDlib_250205a.csv'))
# %% import gsea without doing above, and format for PMF
gsealist = pd.read_csv(os.path.join(filepath3,'gsea_FDR_NES_redux_HBDlib_250205a.csv'), index_col = 0)

# create df with counts, proba and nes
pmf = pd.DataFrame({'count': gsealist['Term.1'].value_counts(),
                    'nes': gsealist.groupby(['Term.1'])['NES'].mean()})

pmf['proba'] = pmf['count']/pmf['count'].sum()

pmf.sort_values(['proba'], inplace = True, ascending = False)

# plot pmf for gsea terms
fig, ax = plt.subplots(figsize =[6,5])
sns.barplot(x = pmf.index, y = pmf['proba'], hue = pmf['nes'], palette = 'RdBu_r')
ax.tick_params(axis='x', labelsize=15, rotation=90)  # Adjust font size and rotation
ax.tick_params(axis='y', labelsize=15)  # Adjust font size
ax.set_xlabel('')  # This removes the x-axis label
ax.set_ylabel('')
plt.savefig(os.path.join(fig_out, 'PMF_GSEA_PROTACs_all.pdf'))

# %% Plot the results
norm = plt.Normalize(vmin=gsealist['NES'].min(), vmax=gsealist['NES'].max())
gradient = 'RdBu_r'

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26'

# Create a stripplot (dot plot)

gseaplot = plt.figure(figsize=(30, 10))

sns.stripplot(x="condition", y="Term.1", data=gsealist, hue="NES",  order = cat_order,
              palette=gradient, jitter=False, size=30, dodge=False)
ax = plt.gca()
ax.yaxis.tick_right()
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gradient), ax=ax, location = 'left')
cbar.set_label('Enrichment')
plt.title("GSEA Cell Compartment, Top " + str(top) + ' by FDR', fontsize = '30')
ax.tick_params(axis='x', labelsize=30, rotation=90)  # Hides x-axis labels
ax.tick_params(axis='y', labelsize = 30)  
ax.set_xlabel('')  # This removes the x-axis label
ax.set_ylabel('') 
ax.get_legend().remove()
plt.show()
gseaplot.savefig(os.path.join(fig_out, 'GSEA_ChemicalSeries_1and10uM.pdf'))

# %% extract t-matrices

matrix_0p1 = pd.read_csv(os.path.join(filepath3, 'Cluster_LFCxPval_0p1uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T

matrix_1 = pd.read_csv(os.path.join(filepath3, 'Cluster_LFCxPval_1uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T

matrix_10 = pd.read_csv(os.path.join(filepath3, 'Cluster_LFCxPval_10uM_250305a.csv'), sep = ';', index_col = 0, decimal = ',').T

title = 'GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX'

toplot = dogma.loc[title]['Lead_genes'][1].split(';')

toplot = dogma.loc[title]['Lead_genes'][1].split(';')

# %%

def format_df(df):
    df1 = df.iloc[[5,14, 11, 0]]
    df1.index = [14, 9, 6, 1]
    return df1

matrix_0p1 = format_df(matrix_0p1)
matrix_0p1['conc'] = 0.1
matrix_1 = format_df(matrix_1)
matrix_1['conc'] = 1
matrix_10 = format_df(matrix_10)
matrix_10['conc'] = 10

m = [matrix_0p1, matrix_1, matrix_10]

m = pd.concat(m)

t_matrix = m[toplot]
t_matrix['conc'] = m['conc']

# %%
plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26'

lineplot2 = plt.figure(figsize = (8,10))

plt.errorbar(t_matrix.loc[14]['conc'], t_matrix.loc[14].drop('conc', axis = 1).mean(axis = 1), t_matrix.loc[14].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'orangered', label = 'CRBN-lenalinomide-5N (cluster 14)', linewidth = 2.5)

plt.errorbar(t_matrix.loc[9]['conc'], t_matrix.loc[9].drop('conc', axis = 1).mean(axis =1), t_matrix.loc[9].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'lightgreen', label = 'CRBN-thalidomide-6N (cluster 9)', linewidth = 2.5)

plt.errorbar(t_matrix.loc[6]['conc'], t_matrix.loc[6].drop('conc', axis = 1).mean(axis =1), t_matrix.loc[6].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'aqua', label = 'CRBN-thalidomide-5N (cluster 6)', linewidth = 2.5)

plt.errorbar(t_matrix.loc[1]['conc'], t_matrix.loc[1].drop('conc', axis = 1).mean(axis = 1), t_matrix.loc[1].drop('conc', axis = 1).sem(axis = 1), marker = 'o', capsize = 3, color = 'darkblue', label = 'VHL-amide (cluster 1)', linewidth = 2.5)

#plt.legend(fontsize = 26)

plt.ylabel('Averaged DE (' + str(t_matrix.shape[1] )+ ' Proteins)')

plt.xlabel('Concentration \u03bcM')

plt.title(title)

plt.show()

lineplot2.savefig(os.path.join(fig_out, title.lower() + '_lineplot.pdf'))

# %%

def extract_genes(list):
    unique_genes = set().union(*list)
    return unique_genes

path_genes =list(extract_genes(dogma2['Lead_genes'].str.split(';')))