# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import brier_score_loss
import joblib
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device_gradientboostingmachine import GBDT 

# %% Relative paths & import datasets

path = os.path.dirname(__file__)

path = os.path.join(path, '..', 'data')

out = os.path.join(path, '..', 'figures')

def clean_drug_index(df):

    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)

    df.index = df.index.str.replace(' - Compound', '')

    return df

LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 'Drug_LFCxPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

expression_matrix =  pd.read_csv(os.path.join(path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T

#  DIA-MS input for Machine Learning

matrix = LFC_matrix.copy()

# %% Format predictor and response vars

AZmeta_input= pd.read_csv(os.path.join(path, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)

AZmeta = AZmeta_input.copy()

AZmeta.index = AZmeta.index.str.replace('-','_') # Format index

AZmeta['Gal'] = AZmeta['Gal'].str.replace('>','').astype(float) # Convert IC50 to float

IC50s = pd.DataFrame({'Gal_IC50': AZmeta['Gal']}, dtype = float) # Extract IC50s
 
IC50s = IC50s.loc[matrix.index] # Match to proteome index

IC50snoNA = IC50s[IC50s['Gal_IC50'].isna()==False] # remove NAs

BinaryTox = pd.DataFrame({ # Make outcome var categorical
    'IC50' : (IC50snoNA['Gal_IC50']<10).astype(int)  
    })

# %% Train-test split

X = matrix.loc[BinaryTox.index].copy() 

y = BinaryTox.copy()

Xb_train, Xb_test, yb_train, yb_test = train_test_split(X, BinaryTox, test_size=0.2, random_state=42) # Baseline Split

AZmeta_input.loc[Xb_train.index.str.replace('_','-'), 'GBDT_split_idx'] = 1

AZmeta_input.to_csv(os.path.join(path, 'AZcompound_metadata_clustered_split_240611a.tsv'))
