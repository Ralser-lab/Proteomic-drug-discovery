# %% filepath and dependencies, load datasets and models
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import joblib
import os
import xgboost as xgb
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

path = os.path.dirname(__file__)

inputout = os.path.join(path, '..', 'data')

modelout = os.path.join(path, '..', 'scoring_models')

figout = os.path.join(path, '..','figures')

modelvalues = pd.read_csv(os.path.join(inputout,'Rplot_Figure4.csv'), index_col = 0)

finalmodel = xgb.XGBClassifier()

finalmodel.load_model(os.path.join(modelout, 'final_model_250305a.json'))

calibratedmodel = joblib.load(os.path.join(modelout, 'final_calibrated_model_250305a.pkl'))

# %% Export for table S3

tableS3 = modelvalues.copy()

tableS3.drop(labels=['Actual', 'Predicted'], axis = 1, inplace = True)

tableS3.columns.values[0] = 'Toxic Probability'

#add Series reference
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14183816').values, 'Drug'] = 'Compound 1'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14196658').values, 'Drug'] = 'Compound 2'
tableS3.loc[pd.Series(tableS3.index).str.contains('AZ14197166').values, 'Drug'] = 'Compound 3'

# check analogue edit
tableS3.loc[tableS3['Drug'].isin(['Compound 1', 'Compound 2', 'Compound 3'])]

tableS3.sort_values('Cluster', inplace = True)

tableS3.to_csv(os.path.join(inputout, 'NCB_ProteomeGuidedDiscovery_TableS3_250606a.csv'))

# %% Cluster 15 and control predictions
cluster15 = modelvalues.loc[modelvalues['Cluster']==15].sort_values('Probability')

targets = ['Moxifloxacin','Troglitazone','Ibuprofen','Ambrisentan', 'Ticrynafen']

controls = modelvalues.loc[modelvalues['Drug'].isin(targets)].sort_values('Probability')

controls.index = controls['Drug']

toplot = pd.concat([cluster15, controls])

#  Merge proba on ID not SN number
values = toplot['Probability'].groupby(toplot.index.str.rsplit('_', n = 1).str[0]).mean().sort_index().sort_values()

# %%
# % Plot toxicity probabilities
plt.rcParams['axes.labelsize'] = '20'   
plt.rcParams['axes.titlesize'] = '20' 
plt.rcParams['xtick.labelsize'] = '20'  
plt.rcParams['ytick.labelsize'] = '20' 
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
plt.figure(figsize=(10, 6))
values.plot(kind = 'bar')
plt.savefig(os.path.join(figout, 'cluster15_predictions.pdf'))
plt.show()

values = list(toplot.index.values)

# %% AZ14183816 analogue predictions 

targets_kevin = 'AZ14183816-005, AZ14183816-006, AZ14196658-003, AZ14197166-003, AZ14197166-004, AZ14197166-005'

targets_kevin = targets_kevin.replace('-', '_').replace(' ','').split(',')

targets_kevin = modelvalues.iloc[:,5:].assign(probability = modelvalues['Probability']).loc[targets_kevin]

#  Protein signatures of compounds

targets_kevin['analogue'] = targets_kevin.index.str.split('_').str[0]

targets_kevin_redux = targets_kevin[['NDUFA5', 'CYC1', 'NDUFA4','NDUFB10', 'NDUFA13', 'analogue']]

plt.rcParams['axes.labelsize'] = '20'   
plt.rcParams['axes.titlesize'] = '20' 
plt.rcParams['xtick.labelsize'] = '20'  
plt.rcParams['ytick.labelsize'] = '20' 
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
plt.figure(figsize=(10, 10))
targets_kevin_redux.plot(kind = 'bar')
plt.savefig(os.path.join(figout, 'targets_signatures.pdf'))
plt.show()

# %% PCA on analogues

analog = targets_kevin.groupby('analogue').mean()

analog.copy().drop('probability', axis = 1, inplace=True)

analog.to_csv(os.path.join(inputout, 'analog_dataout.csv'))

def norm_and_plot(dataframe):
    z = (dataframe - dataframe.mean())/dataframe.std()
    z.plot(kind = 'bar')
    plt.show()
    return z

df = norm_and_plot(analog.loc[:,analog.columns!='probability'])

pca = PCA()
pca.fit(df)
pca.components_

loadings = pd.DataFrame(pca.components_.T, columns=[f'PC_{i+1}' for i in range(len(df.index))], index = df.columns)

components = pd.DataFrame(pca.fit_transform(df), columns=[f'PC_{i+1}' for i in range(len(df.index))], index = df.index)

plt.rcParams['axes.labelsize'] = '20'   
plt.rcParams['axes.titlesize'] = '20' 
plt.rcParams['xtick.labelsize'] = '20'  
plt.rcParams['ytick.labelsize'] = '20' 
plt.rcParams['ytick.labelsize'] = '20'
plt.rcParams['legend.fontsize'] = '20'
plt.rcParams['axes.titlesize'] = '20'
sns.scatterplot(x = 'PC_1', y = 'PC_2',data = components, hue = df.index)

loadings_ratio = pca.components_.T * np.sqrt(pca.explained_variance_)

# Add biplot (arrows) for the loadings
for i, feature in enumerate(loadings.index):
    plt.arrow(0, 0, loadings_ratio[i, 0], loadings_ratio[i, 1], color='red', alpha=0.7, head_width=0.05)
    plt.text(loadings_ratio[i, 0] * 2, loadings_ratio[i, 1] * 1.3, feature, color='black', ha='center', va='center')
plt.legend(title='Analogue', fontsize='small', title_fontsize='medium')  # You can change 'small' to 'x-small' or any size
plt.savefig(os.path.join(figout, 'analogues_signature_PCA.pdf'))
plt.show()

# %% PC2 top features
features = loadings.sort_values(by = 'PC_2').index[0:4]

features = ['NDUFA5', 'CYC1', 'NDUFA4', 'NDUFB10', 'NDUFA13']

analog[features].plot(kind = 'bar', legend = None)
plt.savefig(os.path.join(figout, 'analogues_signature_vs_DMSO.pdf'))

