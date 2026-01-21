#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_21_summarystats_wetlab.py
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
- Figure3_Compound1_Galactose_250204a.xlsx               # Glu/Gal viability titration
- Figure3_compound1_constituents_gal_250604a.xlsx        # Galactose IC50: constituents
- Figure3_Compound1_Seahorse_250204a.xlsx                # Seahorse respiration
- Figure3_all_Complex2_250204a.xlsx                      # Complex II activity (series)
- Figure3_Complex1_250219a.xlsx                          # Complex I activity (Figure 3)
- Figure4_Complex1_250207a.xlsx                          # Complex I activity (Figure 5)
- Figure4_HCC_250208a.xlsx                               # Primary HCC titrations
- Figure4_CellPanel_250314a.xlsx                         # Cell line GI50 panel
- Figure4_XenoGraftARdeg_250207a.xlsx                    # Xenograft AR degradation

Outputs
-------
- figures/protacs_21_compound1_glugal.pdf                       
- figures/protacs_21_galactose.pdf                               
- figures/protacs_21_seahorsemax.pdf                            
- figures/protacs_21_complexII_all.pdf                    
- figures/protacs_21_GI50_dotplot.pdf                        
- figures/protacs_21_xenograft_ARdeg.pdf                         

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
import sys
sys.path.append(os.path.join(os.path.dirname(__file__)))
import device_summarystatistics

# %% Set relative paths
path = os.path.dirname(__file__) 
data_path = os.path.join(path, '..', 'data')
fig_path = os.path.join(path, '..', 'figures')

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

# %% Load, analyze, plot and save glucose galactose compound 1 - 2 - 3 series data. 
condition_list, condition_dfs, norm_dfs = ['cmpd1', 'cmpd2', 'cmpd3', 'untreated'], {}, {}

# Read in all measurements for each treatmenet
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
    plt.title(norm.attrs['title'])
    plt.errorbar(normsum.index, normsum["glu_mean"], yerr=normsum["glu_std"], fmt="o-", label="Glucose Mean", capsize=3, linestyle = '--')
    plt.errorbar(normsum.index, normsum["gal_mean"], yerr=normsum["gal_std"], fmt="s-", label="Galactose Mean", capsize=3, linestyle = '--')
    plt.plot(norm.index, normsum.iloc[:,3:6], color = 'orange', alpha = 0.5)
    plt.plot(norm.index, normsum.iloc[:,0:3], color = 'lightblue', alpha = 0.5)
    plt.ylim(-0.1,1.1)
    plt.savefig(os.path.join(fig_path,f"protacs_21_{norm.attrs['title']}_glugal.pdf"))
    plt.show()
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
plt.tick_params(labelbottom=False, labelleft=False) 
plt.savefig(os.path.join(fig_path, 'protacs_21_galactose.pdf'))
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
plt.savefig(os.path.join(fig_path, 'protacs_21_seahorsemax.pdf'))
plt.close()

# %% Complex II data for all (Figure 3 and Figure 5).

# Read in wetlab source data
complex2 = pd.read_excel(os.path.join(data_path, 'Figure3_all_Complex2_250204a.xlsx'), index_col = 0)

# Perform statistical tests   
device_summarystatistics.one_test(complex2.loc['AZ14183816'], 100)
device_summarystatistics.one_test(pd.Series(complex2.loc[['AZ14197166', 'AZ14183816']].values.flatten()), 100)
 
# Plot and save figure
plot_list = ['AZ14183816', 'AZ14196658','AZ14197166','DMSO']
complex2.loc[plot_list].T.mean().plot(kind = 'bar', 
    yerr = complex2.loc[plot_list].T.std())
plt.savefig(os.path.join(fig_path, 'protacs_21_complexII_all.pdf'))
plt.close()

# %% Perform statistical tests on complex I data for Figure 3.

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

# %% Perform statistical tests on complex I data for Figure 5.

# Read in wetlab source data
complex1 = pd.read_excel(os.path.join(data_path, 'Figure4_Complex1_250207a.xlsx'), index_col = 0)
complex1.index = np.log10(complex1.index) # Numerize to log10index
complex1 = complex1/100 # Bring down to decimal percenet
ic50_values = device_summarystatistics.linear_IC50(complex1) # Calculate log10(IC50) down replicates

# t-test on the estimated log10(IC50s)
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[3:6])
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[6:9])
device_summarystatistics.t_test(ic50_values.iloc[0:3], ic50_values.iloc[3:9])
print((10**ic50_values.iloc[0:3]).mean()) # exponentiate to IC50
print((10**ic50_values.iloc[3:9]).mean()) # exponentiate to IC50

# %% Load and analyze primary HCC data (Figure 5).

# Read in and format wetlab source data
HCC = pd.read_excel(os.path.join(data_path, 'Figure4_HCC_250208a.xlsx'), index_col = 0)
HCC = HCC.loc[~HCC.index.isna()]
HCC2 = HCC.assign(conc = np.log10(HCC.index))
# HCC2.conc.replace(-np.inf, -10, inplace = True)
HCC2 = HCC2.iloc[1:HCC2.shape[0]]
HCC2.set_index('conc', inplace = True)

# Extract column wise IC50s using linear model 
HCC_values = device_summarystatistics.linear_IC50(HCC2)

# Clean up data for t-test
ic50_df = pd.DataFrame({'Compound 1': ic50_values[0:3].values,
                            'Compound 2': ic50_values[3:6].values, 
                            'Compound 3': ic50_values[6:9].values}, 
                            index = ['rep1', 'rep2', 'rep3'])

# Get IC50 values unlogged
ic50s = (10**ic50_df).mean()

# %% Load analyze and plot cell panel data (Figure 5).

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

# Plot and save cell panel data with compound series treatment (Figure 5).
plt.figure(figsize=(6, 6))  # Reduce width
sns.scatterplot(
    data=cell_summ, 
    x='Compound', 
    y='Cell', 
    size=-np.log(cell_summ['mean']), 
    hue=-np.log(cell_summ['mean']), 
    palette='Reds',
    sizes=(50, 500), 
    alpha=0.7
)
plt.xticks(rotation=45, ha='right')  # Rotate labels
plt.xticks(np.arange(0, len(cell_summ['Compound'].unique()), step=1))  # Adjust tick frequency
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.margins(x=0.3, y=0.3) 
plt.tight_layout()  # Adjust layout to reduce excess white space
plt.savefig(os.path.join(fig_path, 'protacs_21_GI50_dotplot.pdf'))
# plt.show()
plt.close()

# %% Perform statistical tests on cell-panel data (Figure 5).

# Melt data for testing
panel = panel.assign(
                factor = pd.Index(panel.index).str.replace(
                r'[.\d]', '', regex=True))
panel = panel.loc[panel['factor']!='Enzalutamide']
ARpos = panel[panel.columns[0:2]]
ARneg = panel[panel.columns[2:panel.shape[1]-2]]

# Perform t-test
device_summarystatistics.t_test(pd.Series(ARpos.values.flatten()),
                pd.Series(ARneg.values.flatten()))

# Perform t-test
ARneg_stat = ARneg.assign(
    factor = ['Compound 1', 'Compound 1', 'Compound 1', 'Compound 1','Analogue', 'Analogue', 'Analogue', 'Analogue', 'Analogue', 'Analogue','Analogue','Analogue'])
device_summarystatistics.t_test(pd.Series(ARneg_stat.loc[ARneg_stat['factor']!= 'Analogue'].iloc[:,0:3].values.flatten()),
                pd.Series(ARneg_stat.loc[ARneg_stat['factor']== 'Analogue'].iloc[:,0:3].values.flatten()))

# %% Load, analyze, plot and save mouse xenograft AR degradation data (Figure 5).

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
plt.savefig(os.path.join(fig_path, 'protacs_21_xenograft_ARdeg.pdf'))
plt.close()



