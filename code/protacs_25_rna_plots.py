# %%
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from functools import lru_cache
from dataclasses import dataclass

@lru_cache
def load_dfs(
    dir_path:str,
    data_fp:str='raw_counts_filtered_logCPM.csv', 
    mart_fp:str='mart_GRCh38.p14.txt',
    meta_fp:str='metadata_limma_design.csv', 
    sets_fp:str='h.all.v2025.1.Hs.symbols.gmt', 
)->Tuple:
    """
    Get the inputs.
    """
    data_df = pd.read_csv(
        os.path.join(
            dir_path, 
            '..',
            'preprocessing_rna',
            'data',
            data_fp,
        ),
        index_col = 0
    )

    meta_df = pd.read_csv(
        os.path.join(
            dir_path, 
            '..',
            'preprocessing_rna',
            'data',
            meta_fp,),
        index_col = 0
    )

    mart_df = pd.read_csv(
        os.path.join(
            dir_path, 
            '..',
            'preprocessing_rna',
            'data',
            mart_fp
        ),
        sep = '\t',
        index_col = 0
    )

    gene_sets = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__), 
            '..',
            'preprocessing_rna',
            'src/gsea/genesets',
            sets_fp,
        ),
        sep = '\t',
        index_col = 0
    )
    return data_df, meta_df, mart_df, gene_sets


def target_search(
    sets:pd.DataFrame, 
    hallmark_id:str,
)->list[str]:
    '''
    Docstring for target_search
    
    :param sets: Description
    :type sets: pd.DataFrame
    :param hallmark_id: Description
    :type hallmark_id: str
    :return: Description
    :rtype: list[str]
    '''
    return sets.loc[hallmark_id].explode().dropna().tolist()

def get_target_ensemblids(
    targets:str,
    mart_df: pd.DataFrame,
)->list[str]:
    '''
    Docstring for get_target_ensemblids
    
    :param targets: Description
    :type targets: str
    :param mart_df: Description
    :type mart_df: pd.DataFrame
    :return: Description
    :rtype: list[str]
    '''
    return (
        mart_df
        .loc[mart_df['Gene name']
        .isin(targets)]
        .index
        .to_list()
    )

def z_condition(
    data_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    ensembl_ids: list[str],
) -> pd.DataFrame:
    
    return (
        data_df
        .pipe(lambda df: (df - df.mean(axis=0)) / df.std(axis=0))
        .loc[:, lambda df: df.columns.isin(ensembl_ids)]
        .assign(condition=meta_df["condition"])
        .groupby("condition")
        .mean()
        .agg(["mean", "sem"], axis=1)                              
    )


def plotter(
        cell_line:str,
        hallmark_id:str,
)->plt.Figure:
    subset_agg_df = z_condition_df.loc[z_condition_df.index.str.contains(cell_line) &
                        ~z_condition_df.index.str.contains('Vehicle|Baseline')]

    df_agg = subset_agg_df.copy()

    # df_agg: index = condition, columns = ['mean', 'sem']
    df = df_agg.reset_index()  # makes a 'condition' column

    # extract timepoint (hours) from condition string
    df["time_h"] = df["condition"].str.extract(r"_(\d+)h").astype(int)

    # assign groups you want as separate lines
    def classify(cond: str) -> str:
        if "Enza" in cond and "3uM" in cond:
            return "Enza (3uM)"
        if "Enza" in cond and "300nM" in cond:
            return "Enza (300nM)"
        if "7166" in cond and "3uM" in cond:
            return "7166 (3uM)"
        if "7166" in cond and "300nM" in cond:
            return "7166 (300nM)"
        return "Other"

    df["group"] = df["condition"].apply(classify)

    # keep only the requested groups and sort by time for proper line plotting
    order = ["7166 (300nM)", "Enza (300nM)", "7166 (3uM)", "Enza (3uM)"]
    order = [g for g in order if g in df["group"].unique()]

    df["group"] = pd.Categorical(df["group"], categories=order, ordered=True)
    df = df[df["group"].isin(order)].sort_values(["group", "time_h"])


    # plot: x = timepoint, separate lines by group, error bars = sem
    fig, ax = plt.subplots()

    for grp in order:
        g = df[df["group"] == grp].sort_values("time_h")
        ax.errorbar(
            g["time_h"],
            g["mean"],
            yerr=g["sem"],
            fmt="-o",
            capsize=4,
            label=grp,
        )

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Pathway activity (z-score ± SEM)")
    ax.set_xticks(sorted(df["time_h"].unique()))
    ax.legend(title="Condition")
    plt.title(f'{cell_line}: {hallmark_id}')
    plt.tight_layout()
    return fig

def plotter_conc(
    cell_line: str,
    hallmark_id: str,
) -> plt.Figure:

    subset_agg_df = z_condition_df.loc[
        z_condition_df.index.str.contains(cell_line) &
        ~z_condition_df.index.str.contains("Vehicle|Baseline")
    ]

    df = subset_agg_df.reset_index(names="condition")

    df["time_h"] = df["condition"].str.extract(r"_(\d+)h").astype(int)

    def get_dose(cond: str) -> str:
        if "300nM" in cond: return "300nM"
        if "3uM" in cond:   return "3uM"
        return None

    def get_drug(cond: str) -> str:
        if "Enza" in cond: return "Enza"
        if "7166" in cond: return "7166"
        return None

    df["dose"] = df["condition"].apply(get_dose)
    df["drug"] = df["condition"].apply(get_drug)
    df = df.dropna(subset=["dose", "time_h", "drug"])

    dose_order = [d for d in ["300nM", "3uM"] if d in df["dose"].unique()]
    df["dose"] = pd.Categorical(df["dose"], categories=dose_order, ordered=True)

    drugs = sorted(df["drug"].unique())
    fig, axes = plt.subplots(1, len(drugs), sharey=True, figsize=(5 * len(drugs), 4))
    if len(drugs) == 1:
        axes = [axes]

    for ax, drug in zip(axes, drugs):
        ddf = df[df["drug"] == drug]

        for t, g in ddf.groupby("time_h"):
            g = g.sort_values("dose")
            ax.errorbar(
                g["dose"],
                g["mean"],
                yerr=g["sem"],
                fmt="-o",
                capsize=4,
                label=f"{t} h",
            )

        ax.set_title(drug)
        ax.set_xlabel("Concentration")
        ax.set_xticks(dose_order)

    axes[0].set_ylabel("Pathway activity \n(z-score logCPM ± SEM)")
    axes[-1].legend(title="Time")

    plt.suptitle(f"{cell_line}: {hallmark_id}")
    plt.tight_layout()
    return fig

def plotter_24h_boxscatter(
    cell_line: str,
    hallmark_id: str,
    time_h: int = 24,
) -> plt.Figure:
    # filter to cell line, exclude Vehicle/Baseline
    subset = z_condition_df.loc[
        z_condition_df.index.str.contains(cell_line) &
        ~z_condition_df.index.str.contains("Vehicle|Baseline")
    ].copy()

    df = subset.reset_index(names="condition")

    # extract time
    df["time_h"] = df["condition"].str.extract(r"_(\d+)h").astype(int)
    df = df[df["time_h"] == time_h].copy()

    # classify into your 4 groups
    def classify(cond: str) -> str:
        if "Enza" in cond and "3uM" in cond:
            return "Enza (3uM)"
        if "Enza" in cond and "300nM" in cond:
            return "Enza (300nM)"
        if "7166" in cond and "3uM" in cond:
            return "7166 (3uM)"
        if "7166" in cond and "300nM" in cond:
            return "7166 (300nM)"
        return None

    df["group"] = df["condition"].apply(classify)
    df = df.dropna(subset=["group"])

    order = ["7166 (300nM)", "Enza (300nM)", "7166 (3uM)", "Enza (3uM)"]
    order = [g for g in order if g in df["group"].unique()]
    df["group"] = pd.Categorical(df["group"], categories=order, ordered=True)
    df = df.sort_values("group")

    # ---- plotting ----
    fig, ax = plt.subplots()

    # boxplot (uses df["mean"] as the value per condition at 24h)
    data = [df.loc[df["group"] == g, "mean"].values for g in order]
    ax.boxplot(data, labels=order, showfliers=False)

    # scatter with jitter
    rng = np.random.default_rng(0)
    for i, g in enumerate(order, start=1):
        y = df.loc[df["group"] == g, "mean"].values
        x = i + rng.normal(0, 0.06, size=len(y))
        ax.scatter(x, y, s=30, alpha=0.8)

    ax.set_ylabel("Pathway activity (z-score logCPM)")
    ax.set_xlabel("Condition")
    ax.set_title(f"{cell_line}: {hallmark_id} @ {time_h} h")
    plt.tight_layout()
    return fig


def plotter(data_df:pd.DataFrame, timepoint, cell, hallmark_id):
    boxplot_df = (
        data_df
        .loc[:, lambda df: df.columns.isin(target_ids)]
        .mean(axis=1)
        .to_frame('Z')
        .assign(condition = meta_df['condition'])
        .loc[lambda df: df['condition'].str.contains(f'{timepoint}')]
        .loc[lambda df: df['condition'].str.contains(f'{cell}')]
        .loc[lambda df: ~df['condition'].str.contains('Baseline')]
    )


    fig, ax = plt.subplots(figsize=(4, 6)) 

    ax = boxplot_df.boxplot(
        column="Z",
        by="condition",
        rot=90,
        showfliers=False,
        ax=ax,
        grid=False
    )

    # scatter overlay
    # for i, cond in enumerate(boxplot_df["condition"].unique(), start=1):
    #     y = boxplot_df.loc[boxplot_df["condition"] == cond, "Z"]
    #     x = np.random.normal(i, 0.04, size=len(y))
    #     ax.scatter(x, y, s=20, alpha=0.7)

    plt.suptitle("")
    plt.title(f"{hallmark_id} @ {timepoint} ({cell})")
    plt.ylabel("Mean target expression (logCPM)")
    plt.tight_layout()

# %%

dir_fp = os.path.dirname(__file__)

data_df, meta_df, mart_df, gene_sets = load_dfs(
    dir_path = dir_fp,
    data_fp = 'raw_counts_filtered_TMM_logCPM.csv',
    meta_fp = 'metadata_limma_design.csv',
    mart_fp = 'mart_GRCh38.p14.txt',
    sets_fp = 'h.all.v2025.1.Hs.symbols.gmt',
)

data_df = data_df.loc[~data_df.index.str.contains('7166_3uM')]
 
hallmark_id:str = ['HALLMARK_ANDROGEN_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION']
timepoints:list[str]= ['6h','24h','72h']
cells:list[str] = ['C4_2']

# %%

from protacs_24_rna_gsea import PlotterCfg, list_getter

cfg = PlotterCfg()

df = list_getter(cfg.dir_name)

targets = df.loc[df['filepath'].str.contains('Enza_3uM_24h')].loc[
    lambda df: df['Term'] == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION'
]['Lead_genes'].str.split(';').explode().tolist()

# %%

for h in hallmark_id: 
    targets = target_search(sets=gene_sets,
        hallmark_id=h)

    target_ids = get_target_ensemblids(targets, mart_df)

    z_condition_df = z_condition(data_df, meta_df, target_ids)

    for c in cells: 
        for t in timepoints:
            plotter(data_df,t, c, h)


# %%

# convert to leading edge for above plots, see if it improves the behavior
targets = target_search(sets=gene_sets,
    hallmark_id='HALLMARK_OXIDATIVE_PHOSPHORYLATION')

len(targets)