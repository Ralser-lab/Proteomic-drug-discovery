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

# %% Galactose IC50 for compound 1 & constituents (Figure 3).

# Read in wetlab sourcedata
cpd1_constituents = pd.read_excel(os.path.join(data_path, 'Figure3_compound1_constituents_gal_250604a.xlsx')).T.iloc[1:,:]
cpd1_constituents.columns = ['AR-binder', 'AR-HBD', 'Lenalidomide']
cpd1_constituents = cpd1_constituents.astype('float64')
order_vec = ['Lenalidomide', 'AR-binder', 'AR-HBD']
cpd1_gal = np.log10(cpd1_constituents)

# Perform test
device_summarystatistics.t_test(cpd1_gal['AR-binder'], cpd1_gal['AR-HBD'])
#device_summarystatistics.t_test(pd.concat([cpd1_gal['AR-binder'],cpd1_gal['Lenalidomide']]), cpd1_gal['AR-HBD'])

# Plot and save figure
plt.figure(figsize=(1,3))
sns.barplot(data = cpd1_gal, order = order_vec, capsize = 0.1, palette = ['grey','red','grey']) #make the ticks better
plt.title('Galactose IC50\nCompound 1 & Constituents', fontsize=13)
plt.tick_params(labelbottom=False, labelleft=False)
plt.savefig(os.path.join(fig_path,f'{workflow}_galactose.pdf'))
plt.close()

# %% Seahorse at max resp (Figure 3).

# Read in wetlab source data
seahorse = pd.read_excel(os.path.join(data_path, 'Figure3_Compound1_Seahorse_250204a.xlsx'))

# Perform t-test
device_summarystatistics.t_test(pd.concat([seahorse['Enza '],seahorse['Control']]), seahorse['Compound 1'])

# Plot and save figure
plt.figure(figsize=(1.5,3))
sns.boxplot(seahorse, order = ['Control','Enza ', 'Compound 1'])
sns.swarmplot(seahorse, color = 'black')
plt.title('Seahorse Max Respiration', fontsize=13)
plt.tick_params(axis='x', labelsize=12)
plt.savefig(os.path.join(fig_path, f'{workflow}_seahorsemax.pdf'))
plt.close()

# %% Complex II data for all (Figure 4 and Figure 6).

# Read in wetlab source data
complex2 = pd.read_excel(os.path.join(data_path, 'Figure3_all_Complex2_250204a.xlsx'), index_col = 0)

# Perform statistical tests   
device_summarystatistics.one_test(complex2.loc['AZ14183816'], 100)
device_summarystatistics.one_test(pd.Series(complex2.loc[['AZ14197166', 'AZ14183816']].values.flatten()), 100)
 
# Plot and save figure
plot_list = ['AZ14183816', 'AZ14196658','AZ14197166','DMSO']
complex2.loc[plot_list].T.mean().plot(kind = 'bar',
    yerr = complex2.loc[plot_list].T.std())
plt.title('Complex II Activity', fontsize=13)
plt.tick_params(axis='x', labelsize=12)
plt.savefig(os.path.join(fig_path, f'{workflow}_complexII_all.pdf'))
plt.close()

# %% Perform statistical tests on complex I data for Figure 4.

# Read in wetlab sourcedata
complex1 = pd.read_excel(os.path.join(data_path, 'Figure3_Complex1_250219a.xlsx'), index_col = 0)
complex1.index = np.log10([0,20,5,1,0.1,0.05,0.01,0.005]) # Numerize to log10 index
complex1 = complex1/100 # Bring down to decimal percent
complex1 = complex1.iloc[1:] # Drop DMSO
ic50_values = device_summarystatistics.linear_IC50(complex1) # calculate log10(IC50) down replicates

# t-test on the estimated log10(IC50s)
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[6:9])
print((10**ic50_values.iloc[0:3]).mean()) # exponentiate to IC50
print((10**ic50_values.iloc[6:9]).mean()) # exponentiate to IC50

# %% Perform statistical tests on complex I data for Figure 6.

# Read in wetlab source data
complex1 = pd.read_excel(os.path.join(data_path, 'Figure4_Complex1_250207a.xlsx'), index_col = 0)
complex1.index = np.log10(complex1.index) # Numerize to log10index
complex1 = complex1/100 # Bring down to decimal percenet
ic50_values = device_summarystatistics.linear_IC50(complex1) # Calculate log10(IC50) down replicates

# t-test on the estimated log10(IC50s)
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[3:6])
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[6:9])
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[3:9])

print('Compound 1 IC50s complex I activity [uM]:')
print((10**ic50_values.iloc[0:3]).mean()) # exponentiate to IC50
print('Analogues IC50s complex I activity [uM]:')
print((10**ic50_values.iloc[3:9]).mean()) # exponentiate to IC50

# %% Glucose-galactose curves and IC50s for compounds 1, 2 and 3 (Figure 6)

# Read in all measurements for each treatment
condition_list, condition_dfs, norm_dfs = ['cmpd1', 'cmpd2', 'cmpd3', 'untreated'], {}, {}
for c in condition_list:
    df = read_sheets(os.path.join(data_path, 'FigureED_Glu_gal_260119a.xlsx'),sheet=c,col=0)
    condition_dfs[c] = df

# Get averages of glucose and galactose untreated controls conditions
glu_untreated_mean = (condition_dfs['untreated']
                      .loc[:, condition_dfs['untreated'].columns.str.contains('Glu')]
                      .mean(axis=1)
                      .iloc[0])
gal_untreated_mean = (condition_dfs['untreated']
                      .loc[:, condition_dfs['untreated'].columns.str.contains('Gal')]
                      .mean(axis=1)
                      .iloc[0])

# Normalize measurements to untreated controls
for c in ['cmpd1', 'cmpd2', 'cmpd3']:
    glu = condition_dfs[c].loc[:,condition_dfs[c].columns.str.contains('Glu')]/glu_untreated_mean
    gal = condition_dfs[c].loc[:,condition_dfs[c].columns.str.contains('Gal')]/glu_untreated_mean
    norm_dfs[c] = pd.concat([glu, gal], axis=1)
    norm_dfs[c].loc[0.00] = 1

# Plot glu-gal curves for each treatment condition
def glu_gal_plotter(norm):
    '''Glu-gal curve plot to disk'''
    normsum = norm.assign(
        glu_mean=norm.iloc[:, 0:3].mean(axis=1),
        glu_std=norm.iloc[:, 0:3].std(axis=1),
        gal_mean=norm.iloc[:, 3:6].mean(axis=1),
        gal_std=norm.iloc[:, 3:6].std(axis=1)
    )
    # Plot and save mean trend lines with error bars
    plt.figure(figsize=(6, 4))
    plt.title(norm.attrs['title'] + ' glu/gal assay', fontsize=13)
    plt.errorbar(normsum.index, normsum["glu_mean"], yerr=normsum["glu_std"], fmt="o-", label="Glucose Mean", capsize=3, linestyle = '--')
    plt.errorbar(normsum.index, normsum["gal_mean"], yerr=normsum["gal_std"], fmt="s-", label="Galactose Mean", capsize=3, linestyle = '--')
    plt.plot(norm.index, normsum.iloc[:,3:6], color = 'orange', alpha = 0.5)
    plt.plot(norm.index, normsum.iloc[:,0:3], color = 'lightblue', alpha = 0.5)
    plt.ylim(-0.1,1.1)
    plt.tick_params(axis='x', labelsize=12)
    plt.savefig(os.path.join(fig_path,f"{workflow}_{norm.attrs['title']}_glugal.pdf"))
    # plt.show()
    plt.close()

for c in norm_dfs:
    glu_gal_plotter(norm_dfs[c])

# Perform statistical test
def IC50_test(df):
    '''
    Performs t-test on IC50s from Glu vs Gal condition 
    estimated from the replicates (n=3).
    
    '''
    real = df.drop([0.00])# Remove normalizing measurement
    print(f"\nGlu:Gal dose-response stats for: {real.attrs['title']}")
    real.index = np.log10(real.index) # Log transform drug concentration
    # Extract replicate (column)-wise IC50s using linear model, clean up for t-test
    ic50_values = device_summarystatistics.linear_IC50(real)
    ic50_df = pd.DataFrame({'Glucose': ic50_values[0:3].values,
                                'Galactose': ic50_values[3:6].values}, 
                                index = ['rep1', 'rep2', 'rep3'])
    # Perform t-test (welch, independent)
    device_summarystatistics.t_test(ic50_df['Glucose'], ic50_df['Galactose'])
    # Unlog IC50
    print('IC50s')
    print(10**(ic50_df.mean()))

for c in norm_dfs:
    IC50_test(norm_dfs[c])

# %% Load and analyze primary HCC data (Figure 6).

# Read in and format wetlab source data
hcc_df = pd.read_excel(os.path.join(data_path, 'Figure4_HCC_250208a.xlsx'), index_col = 0)
hcc_df = hcc_df.loc[~hcc_df.index.isna()] # drop NAs
log_df = hcc_df.assign(conc = np.log10(hcc_df.index)) # log transform X
log_df = log_df.replace([np.inf, -np.inf], 0) # recover 0  
log_df.set_index('conc', inplace=True) 
norm_df = log_df/log_df.max() # normalize Y

# Extract column wise IC50s using linear model 
hcc_ic50s = device_summarystatistics.linear_IC50(norm_df)

# Clean up data for t-test
hcc_ic50_df = pd.DataFrame({'Compound 1': hcc_ic50s[0:3].values,
                            'Compound 2': hcc_ic50s[3:6].values, 
                            'Compound 3': hcc_ic50s[6:9].values}, 
                            index = ['rep1', 'rep2', 'rep3'])

# Get IC50 values unlogged
print('IC50s human hepatocarcinoma survival [uM]:')
print((10**hcc_ic50_df).mean())

# %% Load analyze and plot cell panel GI50 data (Figure 6).

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
vmin, vmax = 0, 30
norm = Normalize(vmin=vmin, vmax=vmax)
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
cbar_ticks = [t for t in [0, 5, 10, 15, 20, 25, 30] if vmin <= t <= vmax]
cbar_tick_labels = [f'{t:.2g}' for t in cbar_ticks]
cbar = fig.colorbar(sm, ax=ax, label='Geomean GI50 (µM)', fraction=0.046, pad=0.04,
                    ticks=cbar_ticks)
cbar.ax.set_yticklabels(cbar_tick_labels)

ax.set_title('Cell Panel GI50', fontsize=13)
ax.tick_params(axis='x', labelsize=12)
plt.margins(x=0.3, y=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_path, f'{workflow}_GI50_dotplot.pdf'))
plt.show()
plt.close()

# %% Load, analyze, plot and save mouse xenograft AR degradation data (Figure 6).

# Read in wetlab source data
degron = pd.read_excel(
                os.path.join(data_path, 'Figure4_XenoGraftARdeg_250207a.xlsx'))

# Perform t-tests
device_summarystatistics.t_test((np.log(degron['Vehicle'])), (np.log(degron['Compound 3'])))
device_summarystatistics.one_test((np.log(degron['Compound 3'])), popmean = np.log(100))

# Plot and save AR degradation in xenograft figure
plt.figure(figsize = [1,3])
sns.boxplot(degron, showfliers = False)
sns.swarmplot(degron, color = 'black')
plt.yscale('log')
plt.title('Xenograft AR Degradation', fontsize=13)
plt.tick_params(axis='x', labelsize=12)
plt.savefig(os.path.join(fig_path, f'{workflow}_xenograft_ARdeg.pdf'))
plt.close()