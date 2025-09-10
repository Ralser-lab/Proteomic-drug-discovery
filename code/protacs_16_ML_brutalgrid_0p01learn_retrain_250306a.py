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

LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 'Drug_LFCxadjPval_250305a.csv'),
                     delimiter = ';', decimal = ',', index_col=0, header = 0).T)

expression_matrix =  pd.read_csv(os.path.join(path,'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'),
                                 delimiter = ',', decimal = '.', index_col=0, header = 0).T

#  DIA-MS input for Machine Learning

matrix = LFC_matrix.copy()

# %% Format predictor and response vars

AZmeta= pd.read_csv(os.path.join(path, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)

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

#  Isolate Enriched Features from split

top5_reactome_enrich_split = pd.read_csv(os.path.join(path, 'top5_FDR_reactome.csv'))

pathway_split = list(top5_reactome_enrich_split.iloc[:,0].value_counts().index)

def extract_genes(df, pathways, matrix, boundary):

    all_genes = set()

    for pathway in pathways:

        path_genes = df[df['term description'] == pathway]['matching proteins in your input (labels)'].values[0].split(',')

        all_genes.update(set(path_genes))

    mask = matrix.mean(axis=0) > boundary

    strong_genes = matrix.columns[mask]

    selected_genes = strong_genes.intersection(all_genes)

    print(len(selected_genes))

    pd.Series(selected_genes).to_csv(os.path.join(path, 'enriched_proteins.csv'), index=0)
    
    return list(selected_genes)

path_genes = extract_genes(top5_reactome_enrich_split, pathway_split, expression_matrix, 11.8)

X_train, X_test, y_train, y_test = train_test_split(X[path_genes], y, test_size=0.2, random_state=42) # Enriched 

# %% first pass GBDT 

dpath = os.path.dirname(__file__)

params = {
                'learning_rate': [0.01],
                'n_estimators': [75,100],
                'max_depth': [1,2,3],
                'min_child_weight': [4,5,6],
                'subsample': [0.8],
                'colsample_bytree': [0.8],
                'gamma': [4,5],
                'alpha': [1,2,3],
                'lambda': [1,2,3]
}

round1 = GBDT(path)

round1.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)

round1.gbdt_gridcv(params, X_train, y_train)

top_features = round1.get_model_features(plot = True, n = X_train.shape[0]//10) # Get top features

# %% second pass GBDT

round2 = GBDT(dpath)

round2.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)

round2.gbdt_gridcv(params, X_train[top_features], y_train)

round2.gbdt_evaluate(X_test[top_features], y_test, round2.best_model)

round2.gbdt_classify(X_test[top_features],y_test, round2.best_model)

round2.gbdt_SHAP(['PRKAR2B','CYC1','NDUFA5','NDUFA4','RPL4','PAFAH1B1','RPL35'])

# %% model calibration

round2.gbdt_calibrate()

round2.gbdt_evaluate(X_test[top_features], y_test, round2.calibrated_model)

round2.gbdt_classify(X_test[top_features], y_test, round2.calibrated_model)

# %% Brier score before and after calibration 

y_pred = round2.calibrated_model.predict(X_test[top_features])

print(brier_score_loss(y_test, y_pred, pos_label = 2))

y_pred = round2.best_model.predict(X_test[top_features])

print(brier_score_loss(y_test, y_pred, pos_label = 2))

# %% export all predictions

round2.gbdt_classify(X[top_features], y['IC50'], round2.calibrated_model)

round2.export_all_predictions(AZmeta, round2.calibrated_model)

# %% export GBDT models

model_out = os.path.join(path, '..' ,'scoring_models')

round2.best_model.get_booster().save_model(os.path.join(model_out, 'final_model_250305a.json')) # final model

joblib.dump(round2.calibrated_model, os.path.join(model_out, 'final_calibrated_model_250305a.pkl')) # calibrated model
