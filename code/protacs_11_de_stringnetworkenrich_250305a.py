# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
import numpy as np
import os

# %% Set up dir

dir_main = os.path.dirname(__file__)

filepath2 = os.path.join(dir_main, '..', 'data/')

filepath3 = os.path.join(dir_main, '..', 'figures')

# %% Import node colors

TXNreact = pd.read_csv(filepath2 + 'Cluster4_enrichment.Component.tsv', sep = '\t')
L5Nreact = pd.read_csv(filepath2 + 'Cluster14_enrichment.Component.tsv', sep = '\t')

# %% Sns striplot

nx_dict = {
    "Nuclear PROTAC": TXNreact,
    "Lenalinomide 5N": L5Nreact,
}

def top5(df):
    df = df.sort_values('false discovery rate', ascending = True)
    df = df.iloc[0:5]
    return df

for n, df in nx_dict.items():
    df = top5(df)
    nx_dict[n] = df

cat_order = ['Nuclear PROTAC', 'Lenalinomide 5N']

dogma = pd.DataFrame()

filter = 0.0001 #FDR cutoff

#Create top GSEA table
for idx, df in nx_dict.items():
    df_subset = df.loc[df['false discovery rate'] < filter]
    df_subset.index = df_subset['term description']
    if df_subset['direction'][0] == 'top':
        df_subset['enrichment score'] = -(df_subset['enrichment score'])
    df_subset['condition'] = idx
    dogma = pd.concat([dogma,df_subset], axis =0)

#Rename pathways to be readable

dogma.loc[dogma['enrichment score'] == dogma['enrichment score'].min(), 'enrichment score'] = dogma['enrichment score'].sort_values()[1]

norm = plt.Normalize(vmin=dogma['enrichment score'].min(), vmax=dogma['enrichment score'].max())
gradient = 'RdBu_r'

custom_palette = ['blue','red'] 

plt.rcParams['axes.labelsize'] = '30'   

plt.rcParams['axes.titlesize'] = '30' 

plt.rcParams['xtick.labelsize'] = '26'  

plt.rcParams['ytick.labelsize'] = '26' 

plt.rcParams['ytick.labelsize'] = '26'

# Create a stripplot (dot plot)

gseaplot = plt.figure(figsize=(10, 10))

sns.stripplot(x="condition", y="term description", data=dogma, hue="enrichment score", order = cat_order,
              palette=gradient, jitter=False, size=30, dodge=False, legend = False)
ax = plt.gca()
ax.yaxis.tick_left()
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=gradient), ax=ax)
cbar.set_label('Enrichment')
plt.title("PPI enrichment, \nTop 5 Pathways (FDR, NES) " + str(filter), fontsize = '30')
ax.tick_params(axis='x', labelsize=18)  
ax.tick_params(axis='y', labelsize = 18)  
ax.set_xlabel('')  
ax.set_ylabel('') 
plt.show()
gseaplot.savefig(os.path.join(filepath3, 'STRINGenrich_PROTACs_OnandOfftarg.pdf'))
