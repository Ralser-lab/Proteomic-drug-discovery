# %%
import os
import pandas as pd

# %%
dirpath = os.path.dirname(__file__)

meta_df = pd.read_csv(
    os.path.join(
        dirpath,
        '../data/sample_metadata_phase2.csv'),
        sep = ','
        )


labels_df = pd.read_csv(
    os.path.join(
        dirpath,
        '../data/251030_Phase2_compound_metadata_withCmax.csv'),
        sep = ','
        )

labels_df
# %%

print(meta_df.columns)
meta_df.head(10)

# %%
def sample_counts(df) -> pd.DataFrame:
    return df.groupby(
        'Compound', 
         dropna = 'False')[[
             'Donor',
             'Dose', 
             'Pool']].describe()

sample_counts(meta_df)

# %% Donor summary

def preturbation_frequencies(df):
    return df.groupby(
    'Dose'
    )['Compound'].describe().sort_values(
    'freq', ascending = False
    ).assign(
    pr_preturb = lambda df:(df['count']/df['count'].sum()*100).astype(int).astype(str)+'%',
    pr_top = lambda df:(df['freq']/df['count']*100).astype(int).astype(str)+'%'
    )

preturbation_frequencies(meta_df)

# %% Donor - Cpd Crosstab

def preturbation_summary(df):
    return pd.crosstab(index = df['Donor'],
    columns=df['Pool']
    )

preturbation_summary(meta_df).style.background_gradient(cmap='viridis')

# %% Concentration - Cpd Crosstab
def replicate_summary(df):
    return df.groupby(
    ['Compound', 'Concentration']
    ).size().unstack(
    fill_value = 0
    ).assign(
    total = lambda df : df.sum(axis = 1)
    )

replicate_summary(meta_df).style.background_gradient(
    cmap = 'viridis'
    )


# %% Dose + Donor - Cpd Crosstab
def replicate_summary(df):
    return df.assign(
        Dose_Donor = df['Concentration'].astype(str) + '_' + df['Donor']
    ).groupby(
    ['Compound', 'Dose_Donor']
    ).size().unstack(
    fill_value = 0
    ).assign(
    total = lambda df : df.sum(axis = 1)
    )

replicate_summary(meta_df).style.background_gradient(
    cmap = 'viridis'
    )


