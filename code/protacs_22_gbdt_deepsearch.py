#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_22_gbdt_deepsearch_train_retrain.py
Description:
    Train and evaluate Gradient Boosted Decision Tree (GBDT) models on 
    proteomic expression data and drug response (IC50) outcomes. The script 
    uses a much deeper hyperpameter search space than protacs_16*.py.

Author: Shaon Basu
Date: 2025-09-29

Inputs
------
- data/Drug_LFCxadjPval_250305a.csv
- data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
- data/AZcompound_metadata_clustered_240611a.tsv
- data/top5_FDR_reactome.csv
- HYPER.json

Outputs
-------
- data/enriched_proteins.csv
- data/predictions.csv 
- figures/*.png
- scoring_models/protacs_22_xgb_first-pass-model.json
- scoring_models/protacs_22_xgb_first-pass-search_space.csv
- scoring_models/protacs_22_xgb_first-pass-final_metrics.csv
- scoring_models/protacs_22_xgb_second-pass-model.json
- scoring_models/protacs_22_xgb_second-pass-search_space.csv
- scoring_models/protacs_22_xgb_second-pass-final_metrics.csv
- scoring_models/protacs_22_xgb_second-pass-calibrated-model.pkl

Requirements
------------
Python >= 3.8  
Dependencies: pandas, scikit-learn, joblib 
Custom: device_gradientboostingmachine 

"""
# %% Import modules
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import brier_score_loss
import joblib
import os
import sys
import json
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device_gradientboostingmachine import GBDT 

# Set relative paths
path = os.path.join(os.path.dirname(__file__), '..', 'data')
model_out = os.path.join(path, '..' ,'scoring_models')
dpath = os.path.dirname(__file__)

# CLI: require HYPER.json
parser = argparse.ArgumentParser(description="GBDT deep hyperparameter search")
parser.add_argument("hyper", help="Path to HYPER.json")
cli_args = parser.parse_args()

# Helper function to format DE matrix index when loading
def clean_drug_index(df):
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

# Load differential expression profiles and proteome from HBD screen
LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 
                                                       'Drug_LFCxadjPval_250305a.csv'
                                                       ), delimiter = ';', decimal = ',', index_col=0, header = 0).T)
expression_matrix =  pd.read_csv(os.path.join(path,
                                              'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'
                                              ),delimiter = ',', decimal = '.', index_col=0, header = 0).T
matrix = LFC_matrix.copy()

# Format predictor and response vars
AZmeta= pd.read_csv(os.path.join(path, 'AZcompound_metadata_clustered_240611a.tsv'),index_col=0)
AZmeta.index = AZmeta.index.str.replace('-','_') # Format index
AZmeta['Gal'] = AZmeta['Gal'].str.replace('>','').astype(float) # Convert IC50 to float
IC50s = pd.DataFrame({'Gal_IC50': AZmeta['Gal']}, dtype = float) # Extract IC50s
IC50s = IC50s.loc[matrix.index] # Match to proteome index
IC50snoNA = IC50s[IC50s['Gal_IC50'].isna()==False] # remove NAs
BinaryTox = pd.DataFrame({ 'IC50' : (IC50snoNA['Gal_IC50']<10).astype(int)}) # Make outcome var categorical

# Format train-test split
X = matrix.loc[BinaryTox.index].copy() 
y = BinaryTox.copy()
Xb_train, Xb_test, yb_train, yb_test = train_test_split(
X, BinaryTox, test_size=0.2, random_state=42) # Baseline Split

# Isolate enriched pathways from split
top5_reactome_enrich_split = pd.read_csv(os.path.join(path, 'top5_FDR_reactome.csv'))
pathway_split = list(top5_reactome_enrich_split.iloc[:,0].value_counts().index)

# Helper function to isolate enriched features
def extract_genes(df, pathways, matrix, boundary):
    """
    For each pathway, the function collects the listed proteins, and saves 
    the final set to 'enriched_proteins.csv'.

    Parameters
    ----------
    df : pandas.DataFrame
    pathways : list of str
    matrix : pandas.DataFrame
    boundary : float

    Returns
    -------
    list of str

    """
    all_genes = set()
    for pathway in pathways:
        path_genes = df[df[
            'term description'] == pathway]['matching proteins in your input (labels)'
                                            ].values[0].split(',')
        all_genes.update(set(path_genes))
    mask = matrix.mean(axis=0) > boundary
    strong_genes = matrix.columns[mask]
    selected_genes = [g for g in matrix.columns if g in strong_genes and g in all_genes]
    pd.Series(selected_genes).to_csv(os.path.join(path, 'enriched_proteins.csv'), index=0)
    return list(selected_genes)

# Isolate enriched features, then split data into training and test sets for machine learning
path_genes = extract_genes(top5_reactome_enrich_split, pathway_split, expression_matrix, 11.8)
X_train, X_test, y_train, y_test = train_test_split(X[path_genes], y, test_size=0.2, random_state=42) 

# %% Gradient boosting workflow
workflow = 'protacs_22'

# Load hyperparameter searchspace from HYPER.json
def _ensure_listify(d):
    """Coerce all values to lists so GridSearchCV accepts them."""
    out = {}
    for k, v in d.items():
        if isinstance(v, list):
            out[k] = v
        else:
            out[k] = [v]
    return out
hyper_path = os.path.abspath(cli_args.hyper)
if not os.path.isfile(hyper_path):
    raise FileNotFoundError(f"HYPER file not found: {hyper_path}")
with open(hyper_path, "r", encoding="utf-8") as f:
    params = json.load(f)
if not isinstance(params, dict):
    raise ValueError("HYPER.json must contain a JSON object (key/value mapping).")
params = _ensure_listify(params)
print(f"[HYPER] Loaded hyperparameters from {hyper_path}")
print(f"[HYPER] Keys: {sorted(params.keys())}")

# First pass GBDT
round1 = GBDT(path, workflow, 'xgb_first-pass')
round1.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
round1.gbdt_gridcv(params, X_train, y_train)
round1.gbdt_evaluate(X_test, y_test, round1.best_model)
round1.gbdt_classify(X_test,y_test, round1.best_model)
top_features = round1.get_model_features(plot = True, n = X_train.shape[0]//10) # Get top features
round1.gbdt_SHAP(top_features)
print(f'Round 1 top features: {top_features}')

# Second pass GBDT
round2 = GBDT(dpath, workflow, 'xgb_second-pass')
round2.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
round2.gbdt_gridcv(params, X_train[top_features], y_train)
round2.gbdt_evaluate(X_test[top_features], y_test, round2.best_model)
round2.gbdt_classify(X_test[top_features],y_test, round2.best_model)
round2.gbdt_SHAP() 

# Model calibration
round2.gbdt_calibrate()
round2.gbdt_evaluate(X_test[top_features], y_test, round2.calibrated_model)
round2.gbdt_classify(X_test[top_features], y_test, round2.calibrated_model)

# Export all predictions
round2.gbdt_classify(X[top_features], y['IC50'], round2.calibrated_model, all = True)
round2.export_all_predictions(AZmeta, round2.calibrated_model) 

# Export GBDT models
round1.best_model.get_booster().save_model(os.path.join( # first pass model
    model_out, f'{round1.dout_prefix}xgb_first-pass-model.json'))
round2.best_model.get_booster().save_model(os.path.join( # second pass model
    model_out, f'{round2.dout_prefix}xgb_second-pass-model.json')) 
joblib.dump(round2.calibrated_model, os.path.join( # calibrated model
    model_out, f'{round2.dout_prefix}xgb_second-pass-calibrated-model.pkl')) 
