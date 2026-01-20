
# %% Import packages
import os
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mygene import MyGeneInfo
import numpy as np
from etl_pipe import etl_pipeline
import seaborn as sns
# # %%
# etl_pipeline()

dirpath = os.path.dirname(__file__)

meta_df = pd.read_csv(
    os.path.join(
        dirpath,
        '../data/sample_metadata_phase2.csv'),
        sep = ','
        )

data_df = pd.read_csv(
    os.path.join(
        dirpath,
        '../data/count_matrix_combat-seq_pool_filtered_TMM_logCPM.csv'),
        sep = ',', index_col = 0
        )

deg_df = pd.read_csv(
    os.path.join(
        dirpath, 
        '../data/deseq2_phase2_compoundwise_corrected_donor_cmpplate_pool_100k_outlier_cutoff_results.csv')
        )

# %% format deg target

cond_df = deg_df.loc[deg_df['Compound']=='caffeine'].assign(
    padjconditioned = lambda df : (df['padj'] < 0.05).astype(int)
)

agg_df = pd.DataFrame(cond_df.groupby('Contrast')['padjconditioned'].sum()).reset_index().assign(
    concentration = lambda df : (df['Contrast'].str.split('_', expand = True)[1]).astype(float),
    proportion = lambda df : df['padjconditioned']#/df['padjconditioned'].sum()*100
).sort_values('concentration')


plt.figure(figsize=(6,5))
plt.bar(agg_df["concentration"].astype(str), agg_df["proportion"])
plt.xlabel("Concentration")
plt.ylabel("Significant Count")

fnt = 14
plt.plot([], [], color='gray', label = 'Contrast =\n~ 0 + Treatment_factor + \nDonor + CmpPlate + Pool')
plt.legend(loc='right', fontsize = fnt)
plt.xlabel('Concentration [uM]', fontsize = fnt)
plt.ylabel('Share of All Caffeine Data (%)', fontsize = fnt)
plt.yticks(fontsize = fnt)
plt.xticks(fontsize = fnt)
plt.title('Caffeine (Significant DEG Counts)', fontsize = fnt)
plt.tight_layout()
plt.show()

cond_df.head()

# %%

#correlation matrix

targets = ['caffeine', 'clotrimazole']

for comp in targets:

    meta_2 = meta_df.set_index('sample_id')

    # Get significant genes
    subset_df = data_df#.loc[:, cond_df.loc[cond_df['padjconditioned'] == 1]['ID']]

    # Subset *only this compound's* samples
    subset_meta = meta_2.loc[meta_2['Compound'] == comp].assign(
        tag=lambda df: df['Compound'] + '_' + df['Donor'] + '_' +
                       df['Concentration'].astype(str) + 'uM'
    )

    subset_df = subset_df.loc[subset_meta.index]

    cleaned_df = subset_df.assign(tag=subset_meta['tag']).set_index('tag', drop=True)

    z_df = (cleaned_df - cleaned_df.mean()) / cleaned_df.std()
    sample_corr = z_df.T.corr()

    g = sns.clustermap(
        sample_corr,
        cmap='RdBu_r',
        figsize=(8, 8),
        row_cluster=True,
        col_cluster=True,
        yticklabels=True,
        xticklabels=True,
        vmin=-1,
        vmax=1
    )

    g.fig.suptitle(comp)  # label which compound this is
    plt.show()
    plt.close(g.fig)


# %%
def format_count_target(data_df:pd.DataFrame, meta_df:pd.DataFrame, target:str)->pd.DataFrame:
    '''Get target in Matrix'''
    set_meta = meta_df[meta_df['Compound']==target]
    set_data = data_df.loc[set_meta['sample_id']]

    # Transform count df
    set_data_sum = pd.DataFrame(
        set_data.var(axis=1), 
        columns=['SumCount']
    ).assign(
        PrSumCount=lambda df: df['SumCount']/ df['SumCount'].sum()*100,
    )

    # Load with meta df for Plotting
    combined_data = pd.merge(set_meta, set_data_sum, left_on='sample_id', right_index = True).assign(
        DoseCat = lambda df : df['Concentration'].astype('category')
    )

    return combined_data

combined_data_caff = format_count_target(subset_df, meta_df, target='caffeine')

# %%

def boxplot_conc(combined_data: pd.DataFrame, 
                 x_type: str='Caffeine (Salmon Aligner Gene Counts)',
                 y_targ: str='PrSumCount'):
    '''Box Plot''' 
    data = combined_data
    doses = sorted(data['Concentration'].unique())
    groups = [data.loc[data['Concentration'] == d, y_targ] for d in doses]

    plt.figure(figsize=(6,5))
    plt.boxplot(groups, labels=doses)

    # Scatter points (donor replicates)
    for i, d in enumerate(doses):
        y = data.loc[data['Concentration'] == d, y_targ]
        x = np.random.normal(i+1, 0.08, size=len(y))
        plt.scatter(x, y, alpha=0.6)

    plt.plot([], [], 'o', color='gray', label='Replicates = Donor')
    plt.legend(loc='upper left', fontsize = fnt)

    plt.xlabel('Concentration [uM]', fontsize = fnt)
    plt.ylabel('Share of All Caffeine Data (%)', fontsize = fnt)
    plt.yticks(fontsize = fnt)
    plt.xticks(fontsize = fnt)
    plt.title(x_type, fontsize = fnt)
    plt.tight_layout()
    plt.show()


boxplot_conc(combined_data_caff)







