#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: protacs_16_gbdt_train_retrain.py
Description:
    Train and evaluate Gradient Boosted Decision Tree (GBDT) models on 
    proteomic expression data and drug response (IC50) outcomes. The script 
    performs a two-stage GBDT workflow with model calibration. Outputs include 
    trained models, feature importance, classification metrics, and prediction tables.

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
from device_gradientboostingmachine import GBDT 
from device_supportfunctions import GBDTUtils
from device_supportfunctions import load_gbdt_inputs
from sklearn.model_selection import train_test_split
from sklearn.metrics import brier_score_loss
import os, joblib

# %% Gradient boosting workflow
def main():

    workflow = 'protacs_17'

    params = {
    'learning_rate': [0.01],
    'n_estimators': [75,100],
    'max_depth': [1,2,3],
    'min_child_weight': [4,5,6],
    'subsample': [0.8],
    'colsample_bytree': [0.8],
    'gamma': [4,5],
    'reg_alpha': [1,2,3],
    'reg_lambda': [1,2,3]
    }

    input = load_gbdt_inputs()
    X, y  = input.X.copy(), input.y.copy()
    Xb_train, Xb_test, yb_train, yb_test = train_test_split(X, y, test_size=0.2, random_state=42) 
    X_train, X_test, y_train, y_test = train_test_split(X[input.path_genes], y, test_size=0.2, random_state=42) 

    # First pass GBDT
    round1 = GBDT(input.path, workflow, 'xgb_first-pass')
    round1.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
    round1.gbdt_gridcv(params, X_train, y_train, score = 'average_precision')
    round1.gbdt_evaluate(X_test, y_test, round1.best_model)
    round1.gbdt_classify(X_test,y_test, round1.best_model) 
    top_features = round1.get_cv_shap_features(plot = True, n = X_train.shape[0]//10) # Get top features
    round1.gbdt_SHAP()

    print('Round1 top features:')
    print(top_features)

    # Second pass GBDT
    round2 = GBDT(input.dpath, workflow, 'xgb_second-pass')
    round2.gbdt_baseline(Xb_train, Xb_test, yb_train, yb_test)
    round2.gbdt_gridcv(params, X_train[top_features], y_train, score = 'average_precision')
    round2.gbdt_evaluate(X_test[top_features], y_test, round2.best_model)
    round2.gbdt_classify(X_test[top_features],y_test, round2.best_model) 
    round2.gbdt_SHAP()

    # Evaluate rounds 
    GBDTUtils.cv_perform(round1, round2, workflow)
    GBDTUtils.plot_pr_two_models(round1.best_model, round2.best_model, X_test, X_test[top_features], y_test, workflow)    

    # Model calibration
    round2.gbdt_calibrate()
    y_pred = round2.best_model.predict(X_test[top_features])
    print(brier_score_loss(y_test, y_pred, pos_label = 2))
    y_pred = round2.calibrated_model.predict(X_test[top_features])
    print(brier_score_loss(y_test, y_pred, pos_label = 2))

    # Export all predictions
    round2.gbdt_classify(X[top_features], y['IC50'], round2.calibrated_model, all = True) # all to True
    round2.export_all_predictions(X[top_features], y['IC50'], input.AZmeta, round2.calibrated_model, workflow) 

    # Export GBDT models
    round1.best_model.get_booster().save_model(os.path.join( # first pass model
        input.model_out, f'{round1.dout_prefix}xgb_first-pass-model.json'))
    round2.best_model.get_booster().save_model(os.path.join( # second pass model
        input.model_out, f'{round2.dout_prefix}xgb_second-pass-model.json')) 
    joblib.dump(round2.calibrated_model, os.path.join( # calibrated model
        input.model_out, f'{round2.dout_prefix}xgb_second-pass-calibrated-model.pkl')) 

if __name__ == "__main__":
    GBDT.configure_font()
    main()


