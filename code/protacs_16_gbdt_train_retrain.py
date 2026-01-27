#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_16_gbdt_train_retrain.py
Description:
    Train and evaluate Gradient Boosted Decision Tree (GBDT) models on 
    proteomic expression data and drug response (IC50) outcomes. The script 
    integrates differential expression (LFC) and Reactome enrichment results, 
    selects enriched pathway features, and performs a two-pass GBDT workflow 
    with model calibration. Outputs include trained models, feature 
    importance, classification metrics, and prediction tables.

Author: Shaon Basu
Date: 2025-09-29

Inputs
------
- data/Drug_LFCxadjPval_250305a.csv
- data/SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv
- data/AZcompound_metadata_clustered_240611a.tsv
- data/top5_FDR_reactome.csv

Outputs
-------
- data/enriched_proteins.csv
- data/predictions.csv 
- figures/*.png
- scoring_models/protacs_16_xgb_first-pass-model.json
- scoring_models/protacs_16_xgb_first-pass-search_space.csv
- scoring_models/protacs_16_xgb_first-pass-final_metrics.csv
- scoring_models/protacs_16_xgb_second-pass-model.json
- scoring_models/protacs_16_xgb_second-pass-search_space.csv
- scoring_models/protacs_16_xgb_second-pass-final_metrics.csv
- scoring_models/protacs_16_xgb_second-pass-calibrated-model.pkl

Requirements
------------
Python >= 3.8  
Dependencies: pandas, scikit-learn, joblib 
Custom: device_gradientboostingmachine, device_supportfunctions

"""
# %% Import modules
import os, sys, joblib
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device_supportfunctions import (
    get_relative_paths,
    get_inputs, 
    format_response_var,
    extract_genes
)
from device_gradientboostingmachine import GBDT 
from sklearn.model_selection import train_test_split
from sklearn.metrics import brier_score_loss

# Set relative paths
path, model_out, dpath = get_relative_paths()

# Load differential expression profiles and proteome from HBD screen
LFC_matrix, expression_matrix, AZmeta = get_inputs(path)

# Format predictor and response vars  
BinaryTox = format_response_var(AZmeta, LFC_matrix)

# Format train-test split with all features
X = LFC_matrix.loc[BinaryTox.index].copy() 
y = BinaryTox.copy()
Xb_train, Xb_test, yb_train, yb_test = train_test_split(
X, BinaryTox, test_size=0.2, random_state=42) # Baseline Split

# Isolate enriched features from split and write to disk 
top5_reactome_enrich_split = pd.read_csv(os.path.join(path, 'top5_FDR_reactome.csv'))
pathway_split = list(top5_reactome_enrich_split.iloc[:,0].value_counts().index)
path_genes = extract_genes(top5_reactome_enrich_split, pathway_split, expression_matrix, 11.8, path)

# Format train-test split with enriched features
X_train, X_test, y_train, y_test = train_test_split(X[path_genes], y, test_size=0.2, random_state=42) 

# %% Gradient boosting workflow
workflow = 'protacs_16'

# Hyperparameter search space
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

# First pass GBDT
round1 = GBDT(path, workflow, 'xgb_first-pass')
round1.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
round1.gbdt_gridcv(params, X_train, y_train)
round1.gbdt_evaluate(X_test, y_test, round1.best_model)
round1.gbdt_classify(X_test,y_test, round1.best_model)
top_features = round1.get_model_features(plot = True, n = X_train.shape[0]//10) # Get top features
round1.gbdt_SHAP(top_features)

print('Round1 top features:')
print(top_features)

# Second pass GBDT
round2 = GBDT(dpath, workflow, 'xgb_second-pass')
round2.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
round2.gbdt_gridcv(params, X_train[top_features], y_train)
round2.gbdt_evaluate(X_test[top_features], y_test, round2.best_model)
round2.gbdt_classify(X_test[top_features],y_test, round2.best_model)
round2.gbdt_SHAP(['PRKAR2B', 'CYC1', 'NDUFA5', 'NDUFA4', 'RPL4', 'PAFAH1B1']) # Top interactors from shap.summary

# Model calibration
round2.gbdt_calibrate()
round2.gbdt_evaluate(X_test[top_features], y_test, round2.calibrated_model)
round2.gbdt_classify(X_test[top_features], y_test, round2.calibrated_model)

# Brier score before and after calibration 
y_pred = round2.best_model.predict(X_test[top_features])
print(brier_score_loss(y_test, y_pred, pos_label = 2))
y_pred = round2.calibrated_model.predict(X_test[top_features])
print(brier_score_loss(y_test, y_pred, pos_label = 2))

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
