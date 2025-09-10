# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import numpy as np
import os

# %% extract input

base = os.path.dirname(__file__)

outpath = os.path.join(base, '..', 'data')

figures = os.path.join(base, '..', 'figures')

# %%

DE_VHL = pd.read_csv(outpath + '/Cluster_VHL_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0, header = 0)

DE_T6N = pd.read_csv(outpath + '/Cluster_T6N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_T5N = pd.read_csv(outpath + '/Cluster_T5N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_L5N = pd.read_csv(outpath + '/Cluster_L5N_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_URA = pd.read_csv(outpath + '/Cluster_URA_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

DE_TXN = pd.read_csv(outpath + '/Cluster_TXN_Limma_250306a.csv',
                     delimiter = ';', decimal = ',', index_col=0)

# %% extract sign log pvalue

DE_list = {
    "DE_VHL": DE_VHL,
    "DE_T6N": DE_T6N,
    "DE_T5N": DE_T5N,
    "DE_URA": DE_URA,
    "DE_L5N": DE_L5N,
    "DE_TXN": DE_TXN,
}

def process_for_network(df):
    df['sign'] = np.sign(df['logFC']) * -np.log10(df['P.Value'])
    ranked = df['sign'].sort_values(ascending=False)
    return ranked

for name, df in DE_list.items():
    ranked = process_for_network(df)
    path = os.path.join(outpath, name + '.csv')
    ranked.to_csv(path)

# %% Sns striplot

# Import reactomes from split dataset

T5Nreact = pd.read_csv(outpath + '/DEsplit_T5N_reactome.tsv', sep = '\t')
TXNreact = pd.read_csv(outpath + '/DEsplit_TXN_reactome.tsv', sep = '\t')
T6Nreact = pd.read_csv(outpath + '/DEsplit_T6N_reactome.tsv', sep = '\t')
L5Nreact = pd.read_csv(outpath + '/DEsplit_L5N_reactome.tsv', sep = '\t')
VHLreact = pd.read_csv(outpath + '/DEsplit_VHL_reactome.tsv', sep = '\t')
URAreact = pd.read_csv(outpath + '/DEsplit_URA_reactome.tsv', sep = '\t')

nx_dict = {
    "Thalidomide 5N": T5Nreact,
    "Transcription": TXNreact,
    "Dihydrouracyl": URAreact,
    "Lenalinomide 5N": L5Nreact,
    "Thalidomide 6N": T6Nreact,
    "VHL amide": VHLreact
}

def top5(df):
    df = df.sort_values('false discovery rate', ascending = True)
    df = df.iloc[0:10]
    return df

for n, df in nx_dict.items():
    df = top5(df)
    nx_dict[n] = df

cat_order = ['Thalidomide 5N','Lenalinomide 5N','VHL amide','Thalidomide 6N', 'Dihydrouracyl','Transcription']

#cat_order = ['Thalidomide 5N','Lenalinomide 5N','VHL amide','Thalidomide 6N','Transcription']

cat_order.reverse() 

dogma = pd.DataFrame()

filter = 0.01 #FDR cutoff

#Create top GSEA table
for idx, df in nx_dict.items():
    df_subset = df.loc[df['false discovery rate'] < filter]
    df_subset.index = df_subset['term description']
    if df_subset['direction'][0] == 'top':
        df_subset['enrichment score'] = -(df_subset['enrichment score'])
    df_subset['condition'] = idx
    dogma = pd.concat([dogma,df_subset], axis =0)

dogma.sort_values('enrichment score', ascending = True)

dogma['term description'] = pd.Categorical(dogma['term description'], ordered = True)

dogma.to_csv(outpath + '/top5_FDR_reactome.csv')
