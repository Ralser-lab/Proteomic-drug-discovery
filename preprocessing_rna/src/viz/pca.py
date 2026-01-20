# %% Import packages
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from etl_pipe import etl_pipeline

# %%

def kde_plot_zstandard(data_df: pd.DataFrame)->tuple:
    figs = []   

    # remove nan genes

    data_df = data_df.dropna(axis=1)

    # remove cols with only 1 unique val in gene

    data_df = data_df.loc[:, data_df.nunique(dropna=True) > 1]

    # First KDE
    fig1 = plt.figure()
    data_df.sample(n=10, axis=1, random_state=42).plot(kind='hist',
                                                       bins = 100,
                                                       density=True)
    plt.title("Raw Data KDE")
    plt.xlabel
    figs.append(fig1)

    # Z-standardize
    zstandard_df = (data_df - data_df.mean()) / data_df.std()

    # Second KDE
    fig2 = plt.figure()
    zstandard_df.sample(n=10, axis=1, random_state=42).plot(kind='kde')
    plt.title("Z-Standardized Data KDE")
    plt.xlim(-2,10)
    figs.append(fig2)

    return zstandard_df, figs


def pca_comp_print(zstandard_df:pd.DataFrame,
                   meta_df:pd.DataFrame)->tuple:
    pca = PCA(n_components=3)

    components = pca.fit_transform(zstandard_df)

    explained_variance = pca.explained_variance_ratio_

    components_df = pd.DataFrame({'PC1': components[:,0],
    'PC2': components[:,1],
    'PC3': components[:,2]}, index = zstandard_df.index)

    components_df = pd.merge(components_df, meta_df, how = 'left', left_index = True, right_index = True)

    return components_df, explained_variance

def plot_3D_PCA(pca_meta_df: pd.DataFrame, 
                key: str, 
                explained_variance, 
                pc1_ax: list = [0.0, 1],
                pc2_ax: list = [0.0, 1],
                pc3_ax: list = [0.0, 1],
                cmap_numeric: str = 'magma_r',
                cmap_cat: str = 'rainbow_r',
                cbar_offset: float = 0.4,
                alpha: float = 1,
                size: int = 40):   
    """
    Plots two 3D scatter plots of PCA results with optional axis swapping for PC2 and PC3.
    """

    components_df = pca_meta_df.copy()

    # ----- Compute axis limits -----
    x_min, x_max = components_df["PC1"].quantile(pc1_ax)
    y_min, y_max = components_df["PC2"].quantile(pc2_ax)
    z_min, z_max = components_df["PC3"].quantile(pc3_ax)

    pad = 0.00

    def add_pad(vmin, vmax):
        rng = vmax - vmin
        if rng == 0:
            rng = 1.0
        return vmin - pad * rng, vmax + pad * rng

    x_min, x_max = add_pad(x_min, x_max)
    y_min, y_max = add_pad(y_min, y_max)
    z_min, z_max = add_pad(z_min, z_max)

    # Update plot style
    plt.rcParams.update({
        'axes.titlesize': 18,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 18
    })

    col = components_df[key]

    # ----- Boolean key branch -----
    if col.dropna().map(lambda x: isinstance(x, bool)).all():
        colors = {True: 'red', False: 'grey'}

        fig = plt.figure(figsize=(18, 8))

        # Plot 1: PC1 vs PC3 vs PC2
        ax1 = fig.add_subplot(121, projection='3d')
        for value in [True, False]:
            subset = components_df[components_df[key] == value]
            ax1.scatter(
                subset['PC1'], subset['PC3'], subset['PC2'],
                color=colors[value], label=f'{key} = {value}',
                alpha=alpha, s=size            
            )
        ax1.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax1.set_ylabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax1.set_zlabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax1.set_title(f'{key} PCA: PC1 vs PC3 vs PC2')
        ax1.legend(title=key)

        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(z_min, z_max)
        ax1.set_zlim(y_min, y_max)

        # Plot 2: PC1 vs PC2 vs PC3
        ax2 = fig.add_subplot(122, projection='3d')
        for value in [True, False]:
            subset = components_df[components_df[key] == value]
            ax2.scatter(
                subset['PC1'], subset['PC2'], subset['PC3'],
                color=colors[value], label=f'{key} = {value}',
                alpha=alpha, s=size            
            )
        ax2.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax2.set_ylabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax2.set_zlabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax2.set_title(f'{key} PCA: PC1 vs PC2 vs PC3')
        ax2.legend(title=key)

        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_zlim(z_min, z_max)

        plt.tight_layout()
        plt.show()

    # ----- Numeric key branch -----
    elif pd.api.types.is_numeric_dtype(col):
        values = col.astype(float)
        vmin, vmax = values.min(), values.max()
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(cmap_numeric)

        fig = plt.figure(figsize=(18, 8))

        # Plot 1
        ax1 = fig.add_subplot(121, projection='3d')
        sc1 = ax1.scatter(
            components_df['PC1'], components_df['PC3'], components_df['PC2'],
            c=values, cmap=cmap, norm=norm, alpha=alpha, s=size 
        )

        ax1.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax1.set_ylabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax1.set_zlabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax1.set_title(f'{key} PCA: PC1 vs PC3 vs PC2')

        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(z_min, z_max)
        ax1.set_zlim(y_min, y_max)

        # Plot 2
        ax2 = fig.add_subplot(122, projection='3d')
        sc2 = ax2.scatter(
            components_df['PC1'], components_df['PC2'], components_df['PC3'],
            c=values, cmap=cmap, norm=norm, alpha=alpha, s=size  
        )

        ax2.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax2.set_ylabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax2.set_zlabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax2.set_title(f'{key} PCA: PC1 vs PC2 vs PC3')

        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_zlim(z_min, z_max)

        cbar = fig.colorbar(sc2, ax=[ax1, ax2], shrink=0.8)
        cbar.set_label(key)

        pos = cbar.ax.get_position()
        cbar.ax.set_position([pos.x0 + cbar_offset, pos.y0, pos.width, pos.height])

        plt.tight_layout()
        plt.show()

    # ----- Categorical key branch -----
    else:
        categories = components_df[key].dropna().unique()
        cmap = plt.get_cmap(cmap_cat)
        num_categories = len(categories)
        colors = {cat: cmap(i / max(num_categories, 1)) for i, cat in enumerate(categories)}

        fig = plt.figure(figsize=(18, 8))

        # Plot 1
        ax1 = fig.add_subplot(121, projection='3d')
        for category in categories:
            subset = components_df[components_df[key] == category]
            ax1.scatter(
                subset['PC1'], subset['PC3'], subset['PC2'],
                color=colors[category], label=str(category),
                alpha=alpha, s=size                   # ← ADDED
            )

        ax1.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax1.set_ylabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax1.set_zlabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax1.set_title(f'{key} PCA: PC1 vs PC3 vs PC2')
        ax1.legend(title=key, bbox_to_anchor=(1.05, 1), loc='upper left')

        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(z_min, z_max)
        ax1.set_zlim(y_min, y_max)

        # Plot 2
        ax2 = fig.add_subplot(122, projection='3d')
        for category in categories:
            subset = components_df[components_df[key] == category]
            ax2.scatter(
                subset['PC1'], subset['PC2'], subset['PC3'],
                color=colors[category], label=str(category),
                alpha=alpha, s=size                   # ← ADDED
            )

        ax2.set_xlabel(f'PC1 ({explained_variance[0] * 100:.1f}%)')
        ax2.set_ylabel(f'PC2 ({explained_variance[1] * 100:.1f}%)')
        ax2.set_zlabel(f'PC3 ({explained_variance[2] * 100:.1f}%)')
        ax2.set_title(f'{key} PCA: PC1 vs PC2 vs PC3')
        ax2.legend(title=key, bbox_to_anchor=(1.05, 1), loc='upper left')

        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_zlim(z_min, z_max)

        plt.tight_layout()
        plt.show()

#%%
        
counts_df = pd.read_csv('../../data/raw_counts.csv', index_col = 0).T
meta_df = pd.read_csv('../../data/metadata.csv', index_col = 0)

# %%
z_standard_df = kde_plot_zstandard(counts_df)[0]

components_df, explained_variance = pca_comp_print(z_standard_df, meta_df)

# %%

for col in meta_df.columns:
    plot_3D_PCA(components_df,
                col,
                explained_variance)