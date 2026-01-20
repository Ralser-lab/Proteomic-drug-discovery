# %%
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Any
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
from dataclasses import dataclass

@dataclass(frozen=True)
class PlotterCfg:
    dir_name: str = "../preprocessing_rna/gsea_output/h_pval/"
    cells: tuple[str, ...] = ("LNCaP", "C4_2")
    priority_terms = [
        "HALLMARK_ANDROGEN_RESPONSE",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
    metrics: tuple[str, ...] = ("NES", "absNES")

def list_files(dir_path: str):
    """Return full paths for all files in directory (skip subfolders/hidden)."""
    return [
        os.path.join(dir_path, f)
        for f in sorted(os.listdir(dir_path))
        if os.path.isfile(os.path.join(dir_path, f)) and not f.startswith(".")
    ]

def sort_gsea_list(file_path:str, topN:int = 20)->pd.DataFrame:
    """Top N |NES| sort of a single gsea output file."""    
    df = pd.read_csv(file_path)

    ranked = (
        df.assign(
            NESabs=lambda d: d["NES"].abs(),
            FDR_q_val=lambda d: pd.to_numeric(d["FDR q-val"], errors="coerce"),
            Gene_prop=lambda d: d["Gene %"].astype(str)
                                .str.replace("%", "", regex=False)
                                .astype(float),
        )
        .loc[lambda d: d["Gene_prop"] != 100]
        .sort_values("NESabs", ascending=False)
        .head(topN)
    )

    return ranked

def list_getter(dir_name:str)->pd.DataFrame:
    df = pd.DataFrame()

    for i in list_files(dir_name):
        df2 = sort_gsea_list(i, 10)
        df2['filepath'] = i
        df = pd.concat([df, df2])

    return df

def clean_filepath(df:pd.DataFrame, target:str)->pd.DataFrame:

    return df.assign(
            filepath=lambda df: (
                df['filepath']
                .str.replace(r'^.*h_pval/', '', regex=True)
                .str.replace('_gsea_results_h.csv', '', regex=False)
                .str.replace('_0', '', regex=False)
                .str.replace(f'{target}_', '', regex=False)
        )
    )

def split_filepath(df:pd.DataFrame)->pd.DataFrame:

    return df.assign(
        conc=lambda df: df['filepath'].str.extract(r'^[^_]+_([^_]+)_').replace(r'\D+', '', regex=True).astype(int),
        time=lambda df: df['filepath'].str.extract(r'_([^_]*)$').replace(r'\D+', '', regex=True).astype(int),
        drug = lambda df: df['filepath'].str.split('_').str[0]
    )

def gsea_summary_plot(
    df: pd.DataFrame,
    target: str = "C4_2",
    priority_terms:list[str]|None=None,
    metric: str = "NES",
    metric2: str = "FDR q-val",
    eps: float = np.finfo(float).eps
) -> Any:

    # ------------- CLEAN 
    # Subset on cellline.
    subset_df = df.loc[df["filepath"].astype(str).str.contains(target, na=False)]

    # Prepare term ordering via priority terms at front followed by frequency.
    topN_terms = subset_df["Term"].value_counts().index.to_list()
    other_terms = [
        t for t in topN_terms
        if t not in priority_terms
    ]
    term_order = priority_terms + other_terms

    # Prepare plotting df.
    nes_df = (
        subset_df.loc[subset_df["Term"].isin(topN_terms)]
        .pipe(clean_filepath, target)
    )

    # Clean filepath col.
    filepath_df = (
        pd.DataFrame({"filepath": nes_df["filepath"].unique()})
        .pipe(split_filepath)
    )

    # Prepare contrast ordering.
    filepath_order = (
        filepath_df.sort_values(["time", "drug", "conc"], ascending=[False, True, True])["filepath"]
        .tolist()
    )

    # Order the plotting df.
    nes_df = (nes_df
              .assign(
                  filepath = lambda df: pd.Categorical(df['filepath'], categories = filepath_order, ordered=True),
                  Term = lambda df: pd.Categorical(df['Term'], categories = term_order, ordered=True),)
              .sort_values(['filepath','Term'])
    )

    # Color normalizer.
    norm = TwoSlopeNorm(
        vmin=nes_df[metric].min(),
        vcenter=0,
        vmax=nes_df[metric].max()
    )

    # ------------- PLOT
    # Size by |NES|.
    nes = nes_df[metric].to_numpy()
    nes_abs = np.abs(nes)
    sizes = 20 + 80 * (nes_abs - nes_abs.min()) / (nes_abs.ptp() + eps)

    # Alpha by -log10(p or q).
    sig = -np.log10(nes_df[metric2].to_numpy() + eps)
    alphas = (sig - sig.min()) / (sig.ptp() + eps)
    alphas = np.clip(alphas, 0.2, 1.0)

    # Color by NES, with per-point alpha via signficance.
    cmap = cm.get_cmap("RdBu_r")
    rgba = cmap(norm(nes))
    rgba[:, 3] = alphas

    x = nes_df["Term"].cat.codes
    y = nes_df["filepath"].cat.codes

    fig, ax = plt.subplots(figsize=(11, 7))

    ax.scatter(x,y,c=rgba,s=sizes)

    # Force axis order from categories.
    ax.set_xticks(np.arange(len(nes_df["Term"].cat.categories)))
    ax.set_xticklabels(nes_df["Term"].cat.categories, rotation=90)

    ax.set_yticks(np.arange(len(nes_df["filepath"].cat.categories)))
    ax.set_yticklabels(nes_df["filepath"].cat.categories)

    ax.set_title(f"Top terms: {target}")
    ax.set_xlabel("Term")
    ax.set_ylabel("Contrast")

    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(nes)
    fig.colorbar(mappable, ax=ax, label=metric)

    fig.tight_layout()
    plt.show()

def gsea_bar_plot(
        df:pd.DataFrame, 
        pathway:str,
        cell:str,
        metric:str,
):
    df_subset = (
        df
        .loc[lambda d: d["Term"].astype(str).str.contains(pathway, na=False, regex=True)]
        .loc[lambda d: d["filepath"].astype(str).str.contains(cell, na=False)]
        .pipe(clean_filepath, cell)
        .assign(sig=lambda d: -np.log10(d["FDR q-val"] + 1e-18),
                absNES = lambda df: np.abs(df['NES']))
        .assign(alpha=lambda d: (d["sig"] - d["sig"].min()) / (d["sig"].max() - d["sig"].min()))
        .assign(alpha=lambda d: d["alpha"].clip(0.2, 1.0))    
        .pipe(split_filepath)
        .set_index('filepath')
        .sort_values(['time', 'drug','conc'], ascending = [True,False,True])
    )

    fig, ax = plt.subplots(figsize=(4, 6))

    bars = ax.bar(df_subset.index, df_subset[metric])

    for bar, a in zip(bars, df_subset["alpha"]):
        bar.set_alpha(float(a))

    ax.set_title(f"{cell}: {pathway}")
    ax.set_ylabel(metric)
    ax.set_xlabel("Contrast")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

# %%
if __name__ == '__main__':
    
    cfg = PlotterCfg()
    gsea_df = list_getter(cfg.dir_name)

    for c in cfg.cells: 
        gsea_summary_plot(
                df=gsea_df,
                target = c,
                metric = 'NES',
                metric2 = 'FDR q-val',
                priority_terms=cfg.priority_terms,
                eps = np.finfo(float).eps)

    for m in cfg.metrics:
        for c in cfg.cells:
            for p in cfg.priority_terms:
                gsea_bar_plot(
                    df=gsea_df,
                    pathway=p,
                    cell=c,
                    metric=m
                    )

