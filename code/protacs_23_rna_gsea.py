#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_23_rna_gsea_plotter.py
Description:
    Load GSEA output CSVs for multiple contrasts, extract the top enriched terms
    by abs NES per file, and generate summary dot plots for selected cell lines.
    Dot color encodes NES, dot size encodes NES, and dot transparency encodes significance
    via -log10(FDR q-val).

Author: Shaon Basu
Date: 2026-01-21

Inputs
------
- ../preprocessing_rna/gsea_output/h_pval/*.csv

Outputs
-------
- figures/protacs_23_rna_<CELL>_gsea_selection_plot.pdf

Requirements
------------
Python >= 3.8
Dependencies: os, numpy, pandas, matplotlib

"""

# %%
import os 
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from typing import Any
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm, Normalize, LinearSegmentedColormap
from dataclasses import dataclass

HERE = os.path.dirname(__file__)

@dataclass(frozen=True)
class PlotterCfg:
    dir_name: str = "../preprocessing_rna/gsea_output/h_pval/"
    cells: tuple[str, ...] = ("LNCaP", "C4_2")
    priority_terms = ["HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_SPERMATOGENESIS",
                      "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_G2M_CHECKPOINT",
                       "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
    metrics: tuple[str, ...] = ("NES", "absNES")
    fig_path_suffix: str = "figures/protacs_23_rna_"

def list_files(dir_path: str)->str:
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

    dir_path = os.path.abspath(os.path.join(HERE, dir_name))

    for i in list_files(dir_path):
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

def gsea_selection_plot(
    df: pd.DataFrame,
    target: str = "C4_2",
    priority_terms:list[str]|None=None,
    metric: str = "NES",
    metric2: str = "FDR q-val",
    eps: float = np.finfo(float).eps
) -> Any:
        
    # Subset on cell-line.
    subset_df = df.loc[df["filepath"].astype(str).str.contains(target, na=False)]

    subset_df = subset_df.loc[~subset_df["filepath"].astype(str).str.contains('7166_3uM',na=False)]
    subset_df = subset_df.loc[~subset_df["filepath"].astype(str).str.contains('Enza_300nM',na=False)]

    # Prepare term ordering via priority terms at front followed by frequency.
    topN_terms = subset_df["Term"].value_counts().index.to_list()

    term_order = priority_terms 

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
              .loc[lambda df: df['Term'].isin(priority_terms)]
    )

    vmin = nes_df[metric].min()
    vmax = nes_df[metric].max()

    if vmin < 0 < vmax:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        cmap = mpl.colormaps["RdBu_r"]
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
        base = plt.cm.Blues
        cmap = LinearSegmentedColormap.from_list(
            "dark_blues",
            base(np.linspace(1.0, 0.8, 256))  # only dark blues
        )

    # Size by |NES|.
    nes = nes_df[metric].to_numpy()
    nes_abs = np.abs(nes)
    sizes = 20 + 80 * (nes_abs - nes_abs.min()) / (nes_abs.ptp() + eps)

    # Alpha by -log10(p or q).
    sig = -np.log10(nes_df[metric2].to_numpy() + eps)
    alphas = (sig - sig.min()) / (sig.ptp() + eps)
    alphas = np.clip(alphas, 0.2, 1.0)

    # Color by NES, with per-point alpha via signficance.
    rgba = cmap(norm(nes))
    rgba[:, 3] = alphas

    x = nes_df["Term"].cat.codes
    y = nes_df["filepath"].cat.codes

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(x,y,c=rgba,s=sizes)

    # Force axis order from categories.
    ax.set_xticks(np.arange(len(nes_df["Term"].cat.categories)))
    ax.set_xticklabels(nes_df["Term"].cat.categories, rotation=90)

    ax.set_yticks(np.arange(len(nes_df["filepath"].cat.categories)))
    ax.set_yticklabels(nes_df["filepath"].cat.categories)

    #ax.set_title(f"Top terms: {target}")
    ax.set_xlabel("Term")
    ax.set_ylabel("Contrast")

    ax.margins(x=0.1)
    ax.margins(y=0.1)

    mpl.rcParams["pdf.fonttype"] = 42   # embed TrueType fonts (Type 42)
    mpl.rcParams["font.family"] = "Helvetica"  # or another installed font

    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(nes)
    fig.colorbar(mappable, ax=ax, label=metric)
    fig.tight_layout()
    plt.savefig(cfg.fig_path_suffix + f'{target}_gsea_dot_plot.pdf')
    #plt.show()
    plt.close()

if __name__ == '__main__':
    
    cfg = PlotterCfg()
    gsea_df = list_getter(cfg.dir_name)

    for c in cfg.cells: 
        
        gsea_selection_plot(
                df=gsea_df,
                target = c,
                metric = 'NES',
                metric2 = 'FDR q-val',
                priority_terms=cfg.priority_terms,
                eps = np.finfo(float).eps)