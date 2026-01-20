# %%
import pandas as pd

meta_df = (
    pd.read_csv('../../data/metadata.csv')
        .assign(
            condition = lambda df: (
            df['cell_line'].str.replace('C4-2','C4_2') + "_" +
            df['drug'] + "_" +
            df['dose'] + "_" +
            df['time'])
        )
        [['cell_line','drug','dose','time','condition']]
        .drop_duplicates()
)
# %%
ref_df = (
    meta_df
    .loc[meta_df['drug']=='Vehicle']
)

targ_df = (
    meta_df
    .loc[~meta_df.index.isin(ref_df.index)&
         ~(meta_df['drug']=='Baseline')]
)

# %%

# build a lookup from (cell_line, time) -> reference condition (Vehicle)
ref_lookup = (
    ref_df
        .set_index(['cell_line', 'time'])['condition']
        .rename('ref_condition')
)

targ_df = (
    targ_df
        .join(ref_lookup, on=['cell_line', 'time'])
        .assign(contrast=lambda d: d['condition'] + " - " + d['ref_condition'])
        .dropna(subset=['ref_condition'])
)

targ_df.to_csv('../../data/contrasts.csv')