# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Any
from matplotlib.colors import TwoSlopeNorm
from config_loader import load_config

# %%
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

def plot_gsea_list(file_path: str, out_dir: str):
    """Plot a dotplot of a single gsea output file."""    
    ranked = sort_gsea_list(file_path)

    if ranked.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    sc = ax.scatter(
        ranked["NES"],          # x = NES
        ranked["Term"],         # y = pathways
        s=300,
        c=-np.log10(ranked["FDR_q_val"] + 1e-6),
        cmap="viridis_r",
        vmin=0,
        vmax=8,
    )

    ax.set_title(os.path.basename(file_path), fontsize=12)
    ax.set_xlabel("NES")
    ax.set_ylabel("Pathway")

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("-log10(FDR q-val)")

    os.makedirs(out_dir, exist_ok=True)
    stem, _ = os.path.splitext(os.path.basename(file_path))
    fig.savefig(os.path.join(out_dir, f"gsea_{stem}.png"), dpi=120)
    plt.close(fig)



def plot_all_gsea_lists(input_dir: str, output_dir: str):
    """Plot multiple dotplots of all gsea output files in a given directory."""   
    for fp in list_files(input_dir):
        plot_gsea_list(fp, output_dir)


if __name__ == "__main__":
    config = load_config()
    plot_all_gsea_lists(
        os.path.join('..',config["output_dir"]), 
        os.path.join('..',config["plot_dir"], 'gsea_by_contrast'))



def gsea_summary_plot(
        dir_name:str, 
        topN:int = 10,
        target:str = 'C4_2',
        metric:str = 'NES',
        metric2:str = 'FDR q-val',
        eps:float = np.finfo(float).eps)->Any:

    df = pd.DataFrame()

    for i in list_files(dir_name):
        df2 = sort_gsea_list(i, topN)
        df2['filepath'] = i
        df = pd.concat([df, df2])

    subset_df = df.loc[df['filepath'].str.contains(target)]

    topN_terms = subset_df['Term'].value_counts().index.to_list()

    nes_df = subset_df.loc[subset_df['Term'].isin(topN_terms)].assign(
        filepath = lambda df: df['filepath']
        .str.replace('../../gsea_output/h_pval/','')
        .str.replace('_gsea_results_h.csv','')
        .str.replace('_0', '')
        .str.replace(f'{target}_', '')
    )

    filepath_df = pd.DataFrame({
        'filepath': nes_df['filepath'].unique()}
        ).assign(
        conc=lambda df: df['filepath'].str.extract(r'^[^_]+_([^_]+)_').replace(r'\D+', '', regex=True).astype(int),
        time=lambda df: df['filepath'].str.extract(r'_([^_]*)$').replace(r'\D+', '', regex=True).astype(int),
        drug = lambda df: df['filepath'].str.split('_').str[0]
    )

    filepath_order = filepath_df.sort_values(['time', 'drug','conc'], ascending = [False,False,True])['filepath'].tolist()

    norm = TwoSlopeNorm(
        vmin=nes_df[metric].min(),
        vcenter=0,
        vmax=nes_df[metric].max()
    )

    nes_df['filepath'] = pd.Categorical(
        nes_df['filepath'],
        categories=filepath_order,
        ordered=True
    )

    nes_df['Term'] = pd.Categorical(
        nes_df['Term'],
        categories=topN_terms,
        ordered=True
    )

    nes_df = nes_df.sort_values(['filepath','Term'])

    sig = -np.log10(nes_df[metric2] + eps)

    sizes = 20 + 80 * (sig - sig.min()) / (sig.max() - sig.min())

    plt.figure(figsize=(11, 7))
    sc = plt.scatter(
        nes_df['Term'],
        nes_df['filepath'],
        c=nes_df[metric],
        s=sizes,
        cmap='RdBu_r',
        norm=norm,
        alpha=0.7
    )
    plt.title(f'Top 10 by |NES|: {target}')
    plt.xlabel('Term')
    plt.xticks(rotation=90)
    plt.ylabel('Contrast')
    plt.colorbar(sc, label=metric)
    plt.tight_layout()
    plt.show()
    return

# %%

gsea_summary_plot('../../gsea_output/h_pval/', 
        topN = 10,
        target = 'C4_2',
        metric = 'NES',
        metric2 = 'FDR q-val',
        eps = np.finfo(float).eps)