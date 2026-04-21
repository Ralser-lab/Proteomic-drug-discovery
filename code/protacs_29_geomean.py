#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_22_summarystats_wetlab.py
Description:
    Analysis and figure generation for wet-lab validation datasets used in Figures 3 and 5. 
    Loads plate/biochemical/cell-panel/xenograft data, computes summary statistics and 
    hypothesis tests, estimates IC50s (log10 scale) with simple linear models, and saves
    publication-ready plots.

Author: Shaon Basu
Date: 2025-09-30

Inputs
------
(data/)
- Figure3_Compound1_Galactose_250204a.xlsx               # Gal viability titration
- FigureED_Glu_gal_260119a                               # Glu/Gal viability titration       
- Figure3_compound1_constituents_gal_250604a.xlsx        # Galactose IC50: constituents
- Figure3_Compound1_Seahorse_250204a.xlsx                # Seahorse respiration
- Figure3_all_Complex2_250204a.xlsx                      # Complex II activity 
- Figure3_Complex1_250219a.xlsx                          # Complex I activity 
- Figure4_Complex1_250207a.xlsx                          # Complex I activity 
- Figure4_HCC_250208a.xlsx                               # Primary HCC titrations
- Figure4_CellPanel_250314a.xlsx                         # Cell line GI50 panel
- Figure4_XenoGraftARdeg_250207a.xlsx  

Outputs
-------
- figures/protacs_22_{cmpd1,cmpd2,cmpd3}_glugal.pdf                       
- figures/protacs_22_galactose.pdf                               
- figures/protacs_22_seahorsemax.pdf                            
- figures/protacs_22_complexII_all.pdf                    
- figures/protacs_22_GI50_dotplot.pdf                        
- figures/protacs_22_xenograft_ARdeg.pdf                         

Requirements
------------
Python >= 3.8
Packages: pandas, numpy, matplotlib, seaborn
Custom: device_summarystatistics 

"""
# %% Import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from matplotlib.colors import Normalize
import device_summarystatistics
from device_supportfunctions import GBDTUtils

# %% Set relative paths
path = os.path.dirname(__file__) 
data_path = os.path.join(path, '..', 'data')
fig_path = os.path.join(path, '..', 'figures')
workflow = 'protacs_22'

GBDTUtils.configure_font()

def add_title(func):
    """Title wrapper for sheet reader. """
    def wrap(*args, **kwargs):
        df = func(*args, **kwargs)
        df.attrs['title'] = kwargs['sheet']
        return df
    return wrap

@add_title
def read_sheets(path, sheet, col):
    """ Read excel sheets. """
    return pd.read_excel(path, sheet_name = sheet, index_col = col)


# %% Load analyze and plot cell panel data (Figure 6).

# Read in and format wetlab sourcedata
panel = pd.read_excel(
                os.path.join(data_path, 'Figure4_CellPanel_250314a.xlsx'),
                index_col = 0).astype(str).T
panel[panel == '>30uM '] = '30000'
panel = panel.astype(float)
panel[panel > 30000] = 30000
panel = panel/1000
panel['Compound'] = ['1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3',
                    'Enz', 'Enz', 'Enz', 'Enz']
panel = panel.loc[panel['Compound']!='Enz']
panel_long=panel.melt(id_vars=['Compound'], var_name =  'Cell', value_name='Value')
AR_pos = ['LNCAP', 'VCAP']
panel_pos = panel_long.loc[panel_long['Cell'].isin(AR_pos)]
panel_neg = panel_long.loc[~panel_long['Cell'].isin(AR_pos)]
cell_summ = panel_long.groupby(['Cell', 'Compound'])['Value'].agg(['mean', 'std']).reset_index()
cell_order = ['LNCAP','VCAP','SU-8686','LU99','MiaPaCa2']
cell_summ['Cell'] = pd.Categorical(cell_summ['Cell'], categories=cell_order, ordered=True)

# %% T-test: LNCAP GI50, Compound 1 vs Compound 3
cell_line = 'SU-8686'

lncap_cpd1 = panel_long.loc[(panel_long['Cell'] == cell_line) & (panel_long['Compound'] == '1'), 'Value']
lncap_cpd3 = panel_long.loc[(panel_long['Cell'] == cell_line) & (panel_long['Compound'] == '3'), 'Value']
print(f'{cell_line} t-test: Compound 1 vs Compound 3')
device_summarystatistics.t_test(lncap_cpd1, lncap_cpd3)

# Boxplot: LNCAP GI50 Compound 1 vs Compound 3
lncap_box = pd.DataFrame({'Cpd 1': lncap_cpd1.values, 'Cpd 3': lncap_cpd3.values})
plt.figure(figsize=(2, 4))
sns.boxplot(lncap_box, showfliers=False)
sns.swarmplot(lncap_box, color='black')
plt.yscale('log')
plt.title(f'{cell_line}', fontsize=13)
plt.ylabel('GI50 (µM)')
plt.tight_layout()
plt.savefig(os.path.join(fig_path, f'{workflow}_LNCAP_GI50_boxplot.pdf'))
plt.show()
plt.close()

# %% Plot and save cell panel data with compound series treatment (Figure 6).
fig, ax = plt.subplots(figsize=(6, 6))
vmin, vmax = round(cell_summ['median'].min(), 2), 30
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sns.scatterplot(
    data=cell_summ,
    x='Compound',
    y='Cell',
    s=300,
    hue=cell_summ['median'],
    palette='Reds_r',
    hue_norm=norm,
    alpha=0.7,
    ax=ax,
    legend=False
)
ax.set_xticks(np.arange(0, len(cell_summ['Compound'].unique()), step=1))
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

# Add a continuous colorbar for mean GI50 values
cmap = plt.get_cmap('Reds_r')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar_ticks = [t for t in [vmin, 5, 10, 15, 20, 25, vmax] if vmin <= t <= vmax]
cbar_tick_labels = [f'{t:.2g}' for t in cbar_ticks]
cbar = fig.colorbar(sm, ax=ax, label='Median GI50 (µM)', fraction=0.046, pad=0.04,
                    ticks=cbar_ticks)
cbar.ax.set_yticklabels(cbar_tick_labels)

ax.set_title('Cell Panel GI50', fontsize=13)
ax.tick_params(axis='x', labelsize=12)
plt.margins(x=0.3, y=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_path, f'{workflow}_GI50_dotplot.pdf'))
plt.show()
plt.close()

# %%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_22_summarystats_wetlab.py
Description:
    Analysis and figure generation for wet-lab validation datasets used in Figures 3 and 5. 
    Loads plate/biochemical/cell-panel/xenograft data, computes summary statistics and 
    hypothesis tests, estimates IC50s (log10 scale) with simple linear models, and saves
    publication-ready plots.

Author: Shaon Basu
Date: 2025-09-30

Inputs
------
(data/)
- Figure3_Compound1_Galactose_250204a.xlsx               # Gal viability titration
- FigureED_Glu_gal_260119a                               # Glu/Gal viability titration       
- Figure3_compound1_constituents_gal_250604a.xlsx        # Galactose IC50: constituents
- Figure3_Compound1_Seahorse_250204a.xlsx                # Seahorse respiration
- Figure3_all_Complex2_250204a.xlsx                      # Complex II activity 
- Figure3_Complex1_250219a.xlsx                          # Complex I activity 
- Figure4_Complex1_250207a.xlsx                          # Complex I activity 
- Figure4_HCC_250208a.xlsx                               # Primary HCC titrations
- Figure4_CellPanel_250314a.xlsx                         # Cell line GI50 panel
- Figure4_XenoGraftARdeg_250207a.xlsx  

Outputs
-------
- figures/protacs_22_{cmpd1,cmpd2,cmpd3}_glugal.pdf                       
- figures/protacs_22_galactose.pdf                               
- figures/protacs_22_seahorsemax.pdf                            
- figures/protacs_22_complexII_all.pdf                    
- figures/protacs_22_GI50_dotplot.pdf                        
- figures/protacs_22_xenograft_ARdeg.pdf                         

Requirements
------------
Python >= 3.8
Packages: pandas, numpy, matplotlib, seaborn
Custom: device_summarystatistics 

"""
# %% Import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from matplotlib.colors import LogNorm
import device_summarystatistics
from device_supportfunctions import GBDTUtils

# %% Set relative paths
path = os.path.dirname(__file__)
data_path = os.path.join(path, '..', 'data')
fig_path = os.path.join(path, '..', 'figures')
workflow = 'protacs_22'

GBDTUtils.configure_font()

def add_title(func):
    """Title wrapper for sheet reader. """
    def wrap(*args, **kwargs):
        df = func(*args, **kwargs)
        df.attrs['title'] = kwargs['sheet']
        return df
    return wrap

@add_title
def read_sheets(path, sheet, col):
    """ Read excel sheets. """
    return pd.read_excel(path, sheet_name = sheet, index_col = col)


# %% Linear

# Read in and format wetlab sourcedata
panel = pd.read_excel(
                os.path.join(data_path, 'Figure4_CellPanel_250314a.xlsx'),
                index_col = 0).astype(str).T
panel[panel == '>30uM '] = '30000'
panel = panel.astype(float)
panel[panel > 30000] = 30000
panel = panel/1000
panel['Compound'] = ['1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3',
                    'Enz', 'Enz', 'Enz', 'Enz']
panel = panel.loc[panel['Compound']!='Enz']
panel_long=panel.melt(id_vars=['Compound'], var_name =  'Cell', value_name='Value')
AR_pos = ['LNCAP', 'VCAP']
panel_pos = panel_long.loc[panel_long['Cell'].isin(AR_pos)]
panel_neg = panel_long.loc[~panel_long['Cell'].isin(AR_pos)]
cell_summ = panel_long.groupby(['Cell', 'Compound'])['Value'].agg(['mean', 'std']).reset_index()
geomean = lambda x: np.exp(np.log(x).mean())
cell_geomean = panel_long.groupby(['Cell', 'Compound'])['Value'].agg(geomean=geomean).reset_index()

key = 'mean'
if key=='geomean':
    cell_summ = cell_geomean.copy()
cell_order = ['LNCAP','VCAP','SU-8686','LU99','MiaPaCa2']
cell_summ['Cell'] = pd.Categorical(cell_summ['Cell'], categories=cell_order, ordered=True)

# Plot and save cell panel data with compound series treatment (Figure 6).
fig, ax = plt.subplots(figsize=(6, 6))
vmin, vmax = 0, 30
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sns.scatterplot(
    data=cell_summ,
    x='Compound',
    y='Cell',
    s=300,
    hue=cell_summ[key],
    palette='Reds_r',
    hue_norm=norm,
    alpha=0.7,
    ax=ax,
    legend=False
)
ax.set_xticks(np.arange(0, len(cell_summ['Compound'].unique()), step=1))
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

# Add a continuous colorbar for mean GI50 values
cmap = plt.get_cmap('Reds_r')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar_ticks = [t for t in [vmin, 5, 10, 15, 20, 25, vmax] if vmin <= t <= vmax]
cbar_tick_labels = [f'{t:.2g}' for t in cbar_ticks]
cbar = fig.colorbar(sm, ax=ax, label=f'{key} GI50 (µM)', fraction=0.046, pad=0.04,
                    ticks=cbar_ticks)
cbar.ax.set_yticklabels(cbar_tick_labels)

ax.set_title('Cell Panel GI50', fontsize=13)
ax.tick_params(axis='x', labelsize=12)
plt.margins(x=0.3, y=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_path, f'{workflow}_GI50_dotplot.pdf'))
plt.show()
plt.close()

# %% Log

# Read in and format wetlab sourcedata
panel = pd.read_excel(
                os.path.join(data_path, 'Figure4_CellPanel_250314a.xlsx'),
                index_col = 0).astype(str).T
panel[panel == '>30uM '] = '30000'
panel = panel.astype(float)
panel[panel > 30000] = 30000
panel = panel/1000
panel['Compound'] = ['1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3',
                    'Enz', 'Enz', 'Enz', 'Enz']
panel = panel.loc[panel['Compound']!='Enz']
panel_long=panel.melt(id_vars=['Compound'], var_name =  'Cell', value_name='Value')
AR_pos = ['LNCAP', 'VCAP']
panel_pos = panel_long.loc[panel_long['Cell'].isin(AR_pos)]
panel_neg = panel_long.loc[~panel_long['Cell'].isin(AR_pos)]
geomean = lambda x: np.exp(np.log(x).mean())
cell_summ = panel_long.groupby(['Cell', 'Compound'])['Value'].agg(geomean=geomean).reset_index()
cell_order = ['LNCAP','VCAP','SU-8686','LU99','MiaPaCa2']
cell_summ['Cell'] = pd.Categorical(cell_summ['Cell'], categories=cell_order, ordered=True)

# Plot and save cell panel data with compound series treatment (Figure 6).
fig, ax = plt.subplots(figsize=(6, 6))
vmin, vmax = 0.1, 30
norm = LogNorm(vmin=vmin, vmax=vmax)
sns.scatterplot(
    data=cell_summ,
    x='Compound',
    y='Cell',
    s=300,
    hue=cell_summ['geomean'],
    palette='Reds_r',
    hue_norm=norm,
    alpha=0.7,
    ax=ax,
    legend=False
)
ax.set_xticks(np.arange(0, len(cell_summ['Compound'].unique()), step=1))
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

# Add a continuous colorbar for geomean GI50 values
cmap = plt.get_cmap('Reds_r')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar_ticks = [t for t in [0.01, 0.1, 1, 10, 30] if vmin <= t <= vmax]
cbar = fig.colorbar(sm, ax=ax, label='Geomean GI50 (µM)', fraction=0.046, pad=0.04)
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels([f'{t:.2g}' for t in cbar_ticks])

ax.set_title('Cell Panel GI50', fontsize=13)
ax.tick_params(axis='x', labelsize=12)
plt.margins(x=0.3, y=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_path, f'{workflow}_GI50_dotplot.pdf'))
plt.show()
plt.close()