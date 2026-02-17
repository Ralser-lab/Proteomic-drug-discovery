#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: device_gradientboostingmachine.py
Description:
    End-to-end Gradient Boosted Decision Trees (GBDT) utilities for binary
    classification using XGBoost. Provides helpers to train a quick baseline,
    run GridSearchCV, evaluate with ROC/PR curves, classify with threshold
    tuning, interpret with SHAP, calibrate probabilities, inspect feature
    importance, generate common plots, and export model outputs for downstream
    analysis (e.g., R plotting).

Author: Shaon Basu
Date: 2025-09-29

Class architecture:
--------
- GBDT
    gbdt_baseline : Train a baseline XGBClassifier.
    gbdt_gridcv : Perform GridSearchCV hyperparameter tuning.
    gbdt_optuna : Perform Bayesian hyperparameter tuning. 
    gbdt_evaluate : Evaluate model with ROC/PR plots.
    gbdt_classify : Classify with adjustable thresholds; report confusion matrices and metrics.
    gbdt_SHAP : Explain predictions with SHAP values and interaction structure.
    gbdt_calibrate : Calibrate predicted probabilities (isotonic).
    get_model_features : Rank and (optionally) plot feature importance.
    plot_* : Visualization utilities (ROC, PR, F1, confusion, summary bars).
    export_all_predictions : Save predictions + SHAP features for all samples.
    export_subset_predictions : Export a focused subset with per-feature barplots.

Requirements
------------
Python >= 3.8  
Dependencies: pandas, numpy, optuna, scikit-learn, xgboost, shap, seaborn, matplotlib

"""
# Import modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import xgboost as xgb
import optuna
import json
from sklearn.model_selection import (GridSearchCV, 
                                     StratifiedKFold, KFold, 
                                     cross_val_score)
from sklearn.metrics import (classification_report, 
                             confusion_matrix, 
                             roc_curve, auc, 
                             f1_score, precision_recall_curve, 
                             ConfusionMatrixDisplay)
from sklearn.calibration import CalibratedClassifierCV, calibration_curve

class GBDT:
    """
    Wrapper around XGBoost classifiers for a binary classification
    workflow. Encapsulates training, tuning, evaluation, calibration,
    interpretation (SHAP), visualization, and export steps, while keeping
    paths and run metadata in one place.

    """
    def __init__(self, dpath, prefix, name):
        """
        Initialize the GBDT class with data/figure paths and identifiers.

        Parameters
        ----------
        dpath : str
            Base directory path for the project. The class will use this to set
            relative paths to `../data` and `../figures`.
        prefix : str
            String prefix to map to figures to generator script (e.g., 'protacs_16*.py').
        name : str
            Name of the object (i.e. 'xgb_firstpass')

        Attributes
        ----------
        dpath : str
            Absolute path to the data directory.
        dout : str
            Absolute path to the figures directory.
        modelout : str
            Absolute path to the models directory. 
        dout_prefix : str
            Prefix for figure filenames.
        name : str
            Name identifier string for the model run.
        baseline_model : xgboost.XGBClassifier or None
            Baseline model trained without tuning.
        best_model : xgboost.XGBClassifier or None
            Best model found after grid search tuning.
        params : dict or None
            Hyperparameter search space for gridCV.
        best_params : dict or None
            Final parameters associated with best model.
        calibrated_model : calibrated xgboost model or None
            Calibrated model using isotonic regression on test set.  
        X_train, y_train : pd.DataFrame, pd.Series
            Training feature matrix and labels (for tuned model).
        X_test, y_test : pd.DataFrame, pd.Series
            Test feature matrix and labels (for tuned model).
        Xb_train, Xb_test, yb_train, yb_test : pd.DataFrame, pd.Series
            Training and test sets for the baseline model.
        proba_test, proba_train : np.ndarray
            Predicted probabilities for test and train sets.
        baseline_proba : np.ndarray
            Predicted probabilities from the baseline model.
        threshold : float
            Decision threshold for classification (default = 0.5).

        """        
        # Relative paths
        self.dpath = os.path.join(dpath, '..', 'data')
        self.dout = os.path.join(dpath, '..', 'figures')
        self.modelout = os.path.join(dpath, '..', 'scoring_models')
        self.dout_prefix = f'{prefix}_'
        self.name = f'{name}_'

        # Late bound
        self.params = None
        self.baseline_model = None
        self.best_params = None
        self.best_model = None
        self.calibrated_model = None
        self.X_train = self.y_train = None
        self.X_test = self.y_test = None
        self.Xb_train = self.Xb_test = self.yb_train = self.yb_test = None
        self.proba_test = self.proba_train = None
        self.baseline_proba = None
        self.threshold = 0.5

    ##############################################################################
    # Helper functions
    ##############################################################################
    
    def out(self, filename, params = False) -> str:
        """
        Build a full output path under the figures directory, prefixed with the
        run's prefix and name. Centralizes figure naming to keep artifacts tidy.

        Parameters
        ----------
        filename : str
            File name (with extension) to append to prefix.
        params : boolean
            If yes, route outpath to /scoring_models

        Returns
        -------
        str
            Full path to the output file in the figures directory.

        """
        if params is False:
            return os.path.join(self.dout, f"{self.dout_prefix}{self.name}{filename}")
        else: 
            return os.path.join(self.modelout, f"{self.dout_prefix}{self.name[:-1]}{filename}")

    ##############################################################################
    # Feature importance 
    ##############################################################################

    def get_model_features(self, plot = True, n = 20) -> pd.Series:
        """
        Extract feature importance (weight, gain, cover) from the fitted best model,
        return the top-n features, and optionally save bar plots for gain and weight.

        Parameters
        ----------
        plot : bool, default True
            If True, save bar plots of top features by Gain and Weight.
        n : int, default 20
            Number of top features to return.

        Returns
        -------
        pandas.Index
            Index of the top-n feature names (columns).

        """
        model = self.best_model
        importances_weight = model.get_booster().get_score(importance_type='weight')
        importances_gain = model.get_booster().get_score(importance_type='gain')
        importances_cover = model.get_booster().get_score(importance_type='cover')

        self.df_importances = pd.DataFrame({'Weight': importances_weight, 
                                            'Gain': importances_gain, 
                                            'Cover': importances_cover}
                                            ).reset_index().rename(columns={'index': 'Feature'})

        self.df_weight_sort = self.df_importances.sort_values(
            ['Weight', 'Gain', 'Cover'], 
            ascending=[False, False, False]).head(n)

        if plot:
            self.plot_top_features(self.df_weight_sort[['Gain', 'Feature']], 'gain')
            self.plot_top_features(self.df_weight_sort[['Weight', 'Feature']], 'weight')

        self.top_n = self.df_importances.sort_values(
            ['Weight', 'Gain', 'Cover'], 
            ascending=[False, False, False]).head(n)['Feature']

        return self.top_n

    ##############################################################################
    # Train / tune / evaluate
    ##############################################################################

    def gbdt_baseline(self, Xb_train, Xb_test, yb_train, yb_test):
        """
        Fit a simple, untuned XGBClassifier as a baseline comparator. Stores the
        model and its test-set predictions/probabilities for side-by-side plots.

        Parameters
        ----------
        Xb_train : pandas.DataFrame
            Baseline training feature matrix.
        Xb_test : pandas.DataFrame
            Baseline test feature matrix.
        yb_train : pandas.Series or numpy.ndarray
            Baseline training labels (binary).
        yb_test : pandas.Series or numpy.ndarray
            Baseline test labels (binary).

        Returns
        -------
        None
            Updates attributes:
            - self.baseline_model
            - self.baseline_predictions
            - self.baseline_proba

        """

        self.Xb_train, self.Xb_test, self.yb_train, self.yb_test = Xb_train, Xb_test, yb_train, yb_test
        baseline_model = xgb.XGBClassifier(random_state=42)
        baseline_model.fit(self.Xb_train, self.yb_train)
        self.baseline_predictions = baseline_model.predict(self.Xb_test)
        self.baseline_proba = baseline_model.predict_proba(self.Xb_test)[:, 1]
        self.baseline_model = baseline_model

    def gbdt_gridcv(self, params, X_train, y_train) -> xgb.XGBClassifier:
        """
        Run a 5-fold GridSearchCV over the provided hyperparameter grid, refit
        the best XGBClassifier on all training data, and store the tuned model.

        Parameters
        ----------
        params : dict
            Grid of hyperparameters for XGBClassifier.
        X_train : pandas.DataFrame
            Training feature matrix.
        y_train : pandas.Series or numpy.ndarray
            Training labels (binary).

        Returns
        -------
        xgboost.XGBClassifier
            The fitted best model found by grid search.

        """
        self.X_train = X_train
        self.y_train = y_train
        self.params = params
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        grid_search = GridSearchCV(
            estimator=xgb.XGBClassifier(),
            param_grid=params,
            cv=kf,
            scoring='accuracy',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)
        best_params = grid_search.best_params_
        best_model = xgb.XGBClassifier(**best_params, random_state=42)
        best_model.fit(X_train, y_train)
        self.best_model = best_model
        self.best_params = best_params
        print(params)
        print(best_params)
        params_df = pd.DataFrame(list(params.items()), columns = ['parameter', 'values'])
        params_df.to_csv(self.out('-search_space.csv', params = True))
        best_params_df = pd.DataFrame(list(best_params.items()), columns = ['parameter', 'value'])
        best_params_df.to_csv(self.out('-final_metrics.csv', params = True))
        return best_model
    
    def gbdt_optuna(self, params, X_train, y_train, n_trials=50, score='accuracy',
                timeout=None, random_state=42, enqueue_params=None) -> xgb.XGBClassifier:
        """
        Run Bayesian search over the provided hyperparameter grid using 5 stratified-splits
        then refit the best XGBClassifier on all training data, and store the tuned model.

        Parameters
        ----------
        params : dict
            Grid of hyperparameters for XGBClassifier.
        X_train : pandas.DataFrame
            Training feature matrix.
        y_train : pandas.Series or numpy.ndarray
            Training labels (binary).

        Returns
        -------
        xgboost.XGBClassifier
            The fitted best model found by grid search.

        """
        self.X_train = X_train
        self.y_train = y_train
        self.params = params

        kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)

        def objective(trial: optuna.Trial) ->float:
            tp = {}
            for name, spec in params.items():
                t = spec["type"]
                if t == "float":
                    tp[name] = trial.suggest_float(
                        name, float(spec["low"]), float(spec["high"]),
                        log=(spec.get("scale") == "log")
                    )
                elif t == "int":
                    tp[name] = trial.suggest_int(
                        name, int(spec["low"]), int(spec["high"]),
                        step=int(spec.get("step", 1))
                    )
                elif t == "categorical":
                    tp[name] = trial.suggest_categorical(name, list(spec["choices"]))
                else:
                    raise ValueError(f"Unknown type for {name}: {t}")

            tp.setdefault("random_state", random_state)
            tp.setdefault("eval_metric", "logloss")
            tp.setdefault("n_jobs", -1)

            model = xgb.XGBClassifier(**tp)
            scores = cross_val_score(model, X_train, y_train, cv=kf, scoring=score, n_jobs=-1)
            return float(np.mean(scores))

        study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(seed=random_state)
        )

        # Enqueue a known-good starting point (optional)
        if enqueue_params is not None:
            # Only pass params that are actually in your search space
            filtered = {k: v for k, v in enqueue_params.items() if k in params}
            study.enqueue_trial(filtered)

        study.optimize(objective, n_trials=n_trials, timeout=timeout, show_progress_bar=False)

        best_params = study.best_trial.params.copy()
        best_params.setdefault("random_state", random_state)
        best_params.setdefault("eval_metric", "logloss")
        best_params["n_jobs"] = -1

        best_model = xgb.XGBClassifier(**best_params)
        best_model.fit(X_train, y_train)

        self.best_model = best_model
        self.best_params = best_params

        pd.DataFrame([(k, json.dumps(v)) for k, v in params.items()], columns=["parameter", "values"])\
        .to_csv(self.out("-search_space.csv", params=True), index=False)
        pd.DataFrame(list(best_params.items()), columns=["parameter", "value"])\
        .to_csv(self.out("-final_metrics.csv", params=True), index=False)

        print("[Optuna] best value:", study.best_value)
        print("[Optuna] best params:", best_params)
        return best_model


    def gbdt_evaluate(self, X_test, y_test, model):
        """
        Evaluate a model on held-out data. Computes predicted probabilities and
        saves ROC and Precision-Recall curves for test (and optionally train) sets,
        including an overlay of the baseline model if available.

        Parameters
        ----------
        X_test : pandas.DataFrame
            Test feature matrix.
        y_test : pandas.Series or numpy.ndarray
            Test labels (binary).
        model : xgboost.XGBClassifier
            Trained model to evaluate.

        Returns
        -------
        None
            Saves ROC and PR plots; stores probabilities for later steps.

        """
        self.X_test = X_test
        self.y_test = y_test
        self.model = model
        self.proba_test = self.model.predict_proba(X_test)[:, 1]
        self.proba_train = self.model.predict_proba(self.X_train)[:, 1]
        self.plot_roc(y_test=y_test, proba=self.proba_test, baseline=self.baseline_proba,
                      y_train=self.y_train, proba2=self.proba_train)
        self.plot_pr(y_test=y_test, proba=self.proba_test, baseline=self.baseline_proba,
                     y_train=self.y_train, proba2=self.proba_train)

    def gbdt_classify(self, X, y, model='best', threshold='auto', all = False):
        """
        Generate binary classifications from predicted probabilities using either
        the stored best model or a provided model. Select a decision boundary
        automatically (max F1) or use a specified threshold, then save confusion
        matrices and a compact classification-report summary plot.

        Parameters
        ----------
        X : pandas.DataFrame
            Feature matrix to classify.
        y : pandas.Series or numpy.ndarray
            True labels for evaluation.
        model : {'best'} or xgboost.XGBClassifier, default 'best'
            Which model to use for prediction. If 'best', uses self.best_model.
        threshold : {'auto', '} or float, default 'auto'
            If 'auto', selects threshold that maximizes F1 on y vs proba.
            If a float in [0, 1], uses that as the decision boundary.
        all : boolean
            If true, use stored decision boundary (F1 score max)

        Returns
        -------
        None
            Saves confusion matrices and classification report plot; updates:
            - self.threshold
            - self.prediction, self.predict_test, self.predict_train, self.predict_baseline
            - self.report

        """
        self.X = X
        self.y = y

        if model == 'best':
            model = self.best_model
            self.proba = model.predict_proba(self.X)[:, 1]
            if threshold == 'auto':
                self.threshold = self.plot_f1_and_pick_threshold(self.y, self.proba, xline=0.50)
            else:
                self.plot_f1_and_pick_threshold(self.y, self.proba, xline=float(threshold))
                self.threshold = float(threshold)
        else:
            self.proba = model.predict_proba(self.X)[:, 1]
            if all is not True:
                if threshold == 'auto':
                    self.threshold = self.plot_f1_and_pick_threshold(self.y, self.proba, xline=0.50)
                else:
                    self.plot_f1_and_pick_threshold(self.y, self.proba, xline=float(threshold))
                    self.threshold = float(threshold)

        self.prediction = self.proba >= self.threshold
        self.predict_test = (self.proba_test >= self.threshold)
        self.predict_train = (self.proba_train >= self.threshold)
        self.predict_baseline = (self.baseline_proba >= self.threshold)

        if all is not True:
            self.plot_confusion(self.y_test, self.proba_test, name='test', boundary=self.threshold)
            self.plot_confusion(self.y_train, self.proba_train, name='train', boundary=self.threshold)
            self.plot_confusion(self.y, self.proba, name='all', boundary=self.threshold)

        print(classification_report(self.y_test, self.predict_test))

        def format_report(df, idx):
            df['accuracy'] = df.loc['accuracy'].mean()
            df = df.iloc[:, [0, 1, 4]].iloc[[1]]
            df.rename(index={'1': idx}, inplace=True)
            df = df * 100
            return df

        final_report = format_report(pd.DataFrame(
            classification_report(self.y_test, self.predict_test, output_dict=True)
        ).T, 'final test')

        baseline_report = format_report(pd.DataFrame(
            classification_report(self.yb_test, self.predict_baseline, output_dict=True)
        ).T, 'baseline test')

        train_report = format_report(pd.DataFrame(
            classification_report(self.y_train, self.predict_train, output_dict=True)
        ).T, 'final train')

        self.report = pd.concat([baseline_report, final_report, train_report])
        self.plot_classification_report(self.report, round(self.threshold * 100))

    def gbdt_SHAP(self, top_interactors=None):
        """
        Compute SHAP values for the best model to quantify feature influence,
        save a SHAP summary plot, derive top features, and visualize SHAP
        interaction structure (overall and within selected 'top_interactors').

        Parameters
        ----------
        top_interactors : list of str or None, default None
            Subset of features to display in the interaction heatmap.
            If None, top interactors from shap.TreeExplainer is used.

        Returns
        -------
        None
            Saves SHAP summary and interaction visualizations to disk; updates:
            - self.df_shap_values
            - self.shap_feature_importance
            - self.top_interactors

        """
        self.explainer = shap.TreeExplainer(self.best_model)
        self.shap_values = self.explainer.shap_values(self.X_test)

        shap.summary_plot(self.shap_values, self.X_test, show=False, plot_size=[5, 10])
        plt.savefig(self.out('shap_explainer.pdf'))
        plt.close()

        self.shap_abs_mean = np.abs(self.shap_values).mean(axis=0)
        self.shap_feature_importance = pd.DataFrame({
            'feature': self.X_test.columns,
            'importance': self.shap_abs_mean
        }).sort_values(by='importance', ascending=False).head(20)

        self.df_shap_values = pd.DataFrame(self.shap_values, index=self.X_test.index, columns=self.X_test.columns)

        interaction_values = self.explainer.shap_interaction_values(self.X_train)
        shap.summary_plot(interaction_values, self.X_train, show=False)
        plt.savefig(self.out('shap_interactions.pdf'))
        plt.close()

        if top_interactors is None:
            top_interactors = self.X_train.columns
        self.top_interactors = top_interactors

        variances = np.std(interaction_values, axis=0)
        variance_df = pd.DataFrame(variances, index=self.X_train.columns, columns=self.X_test.columns)
        filtered_variance_df = variance_df.loc[top_interactors, top_interactors]

        plt.figure(figsize=(10, 8))
        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        sns.clustermap(filtered_variance_df, annot=False, cmap='magma_r', linewidths=.5)
        plt.title('SHAP Interactions in Test set')
        plt.savefig(self.out('shap_interactions_heatmap.pdf'))
        # plt.show()
        plt.close()

    def gbdt_calibrate(self):
        """
        Calibrate class probabilities for both best and baseline models using
        isotonic regression on the test set, then plot and save a calibration
        curve comparing pre- and post-calibration reliability.

        Parameters
        ----------
        None

        Returns
        -------
        None
            Updates attributes:
            - self.calibrated_model
            - self.baseline_predictions
            - self.baseline_proba
            - self.name 
            Saves 'calibration_curve.pdf'.

        """
        calibrated_model = CalibratedClassifierCV(self.best_model, method='isotonic', cv='prefit')
        calibrated_model.fit(self.X_test, self.y_test)

        calibrated_baseline = CalibratedClassifierCV(self.baseline_model, method='isotonic', cv='prefit')
        calibrated_baseline.fit(self.Xb_test, self.yb_test)
        self.baseline_predictions = calibrated_baseline.predict(self.Xb_test)
        self.baseline_proba = calibrated_baseline.predict_proba(self.Xb_test)[:, 1]

        prob_pos = calibrated_model.predict_proba(self.X_test)[:, 1]
        fraction_before, mean_predicted_before = calibration_curve(self.y_test, self.proba_test, n_bins=10)
        fraction_of_positives, mean_predicted_value = calibration_curve(self.y_test, prob_pos, n_bins=10)

        plt.figure(figsize=(10, 10))
        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26 
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        plt.plot(fraction_before, mean_predicted_before, "s-", label="Not Calibrated")
        plt.plot(mean_predicted_value, fraction_of_positives, "s-", label="Calibrated")
        plt.plot([0, 1], [0, 1], "k--", label="Perfectly calibrated")
        plt.xlabel("Mean predicted probability")
        plt.ylabel("Fraction of positives")
        plt.legend()
        plt.savefig(self.out('calibration_curve.pdf'))
        plt.close()
        # plt.show()
        
        self.calibrated_model = calibrated_model
        self.name = self.name[:-1] + '-calibrated_'

    ##############################################################################
    # Export functions
    ##############################################################################

    def export_all_predictions(self, AZmeta, model):
        """
        Export per-sample predictions and probabilities (plus SHAP feature columns)
        merged with metadata (Drug ID, Cluster). Produces a tidy CSV for R
        visualization and stores an in-memory table for further use.

        Parameters
        ----------
        AZmeta : pandas.DataFrame
            Metadata table indexed like X/y with at least ['Drug ID','Dend'].
        model : xgboost.XGBClassifier or CalibratedClassifierCV
            The classifier used to generate probabilities.

        Returns
        -------
        None
            Saves 'Rplot_Figure4.csv' to the data path; updates:
            - self.outcome
            - self.dotplotforR

        """
        self.proba = model.predict_proba(self.X)[:, 1]
        self.prediction = self.proba >= self.threshold
        self.outcome = pd.DataFrame({'Actual': self.y,
                                     'Predicted': self.prediction,
                                     'Probability': self.proba})

        subset = self.outcome.copy()
        subset['Actual'].replace({1: 'Toxic', 0: 'Non-Toxic'}, inplace=True)
        subset['Drug'] = AZmeta.loc[subset.index]['Drug ID']
        subset['Cluster'] = AZmeta.loc[subset.index]['Dend']
        subset.sort_values(['Cluster'])

        clusters = subset.sort_values(['Drug', 'Actual'], ascending=[True, False])
        cluster_shap = self.X.loc[clusters.index][self.shap_feature_importance['feature'].values]
        dotplotforR = pd.merge(clusters, cluster_shap, left_index=True, right_index=True)

        self.dotplotforR = dotplotforR
        dotplotforR.to_csv(os.path.join(self.dpath, 'Rplot_Figure4.csv'))

    def export_subset_predictions(self, AZmeta, targets):
        """
        Export predictions for a subset of 'targets' and generate per-feature
        differential expression bar plots for those targets. Also prepares a
        SHAP-feature-augmented table for R plotting.

        Parameters
        ----------
        AZmeta : pandas.DataFrame
            Metadata table indexed like X/y with at least ['Drug ID','Dend'].
        targets : list of str
            Subset of index labels (e.g., compound IDs) to include.

        Returns
        -------
        None
            Saves individual feature plots and 'Rplot_Figure5.csv'; updates:
            - self.dotplotforR_subset

        """
        subset = self.outcome.copy()
        subset['Actual'].replace({1: 'Toxic', 0: 'Non-Toxic'}, inplace=True)
        subset['Drug'] = AZmeta.loc[subset.index]['Drug ID']
        subset['Cluster'] = AZmeta.loc[subset.index]['Dend']
        subset.sort_values(['Cluster'])

        clusters = subset.sort_values('Drug')
        clusters = clusters.loc[targets]

        self.get_model_features(plot=False)
        gene = list(self.df_importances.sort_values('Weight').iloc[-20:]['Feature'])

        matrix = self.X.copy()
        matrix_subset = matrix[gene]

        def plot_diffexp(column: pd.Series):
            """
            Plot a single feature's sorted values for the selected targets.

            Parameters
            ----------
            column : pandas.Series
                A column (feature) slice across selected targets.

            Returns
            -------
            None
            """
            sorted_col = column.sort_values()
            plt.rcParams['axes.labelsize'] = 35
            plt.rcParams['axes.titlesize'] = 35
            plt.rcParams['xtick.labelsize'] = 20
            plt.rcParams['ytick.labelsize'] = 26
            plt.rcParams['legend.fontsize'] = 28
            plt.figure(figsize=[10, 10])
            plt.title(f'final model weight (F-1): {column.name}\n (AR PROTACs, Thalidomide 5N)')
            plt.ylabel('test set targets vs DMSO')
            sorted_col.plot(kind='bar')
            plt.savefig(self.out(f'{column.name}_LFCxPval.pdf'))
            plt.close()
            # plt.show()

        matrix_subset.loc[targets].apply(plot_diffexp, axis=0)

        cluster_shap = self.X.loc[clusters.index][self.shap_feature_importance['feature'].values]
        dotplotforR = pd.merge(clusters, cluster_shap, left_index=True, right_index=True)

        self.dotplotforR_subset = dotplotforR
        dotplotforR.to_csv(os.path.join(self.dpath, 'Rplot_Figure5.csv'))

    ##############################################################################
    # Plotting functions
    ##############################################################################

    def plot_roc(self, y_test, proba, baseline=None, y_train=None, proba2=None):
        """
        Save ROC curves comparing the evaluated model on test (and optional train)
        data, with an optional overlay of the baseline model evaluated on test.

        Parameters
        ----------
        y_test : pandas.Series or numpy.ndarray
            True labels for the test set.
        proba : numpy.ndarray
            Predicted probabilities on X_test.
        baseline : numpy.ndarray or None, default None
            Baseline model probabilities on the same test labels, if available.
        y_train : pandas.Series or numpy.ndarray or None, default None
            Training labels.
        proba2 : numpy.ndarray or None, default None
            Predicted probabilities on X_train .

        Returns
        -------
        None
            Saves 'roc.pdf' to figures directory.

        """
        fpr, tpr, _ = roc_curve(y_test, proba)
        roc_auc = auc(fpr, tpr)
        if y_train is not None and proba2 is not None:
            fpr_train, tpr_train, _ = roc_curve(y_train, proba2)
            roc_auc_train = auc(fpr_train, tpr_train)

        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        plt.figure(figsize=(10, 10))
        plt.plot(fpr, tpr, color='blue', lw=2, label=f'final (test) \nAUC = {roc_auc:.2f}')
        if y_train is not None and proba2 is not None:
            plt.plot(fpr_train, tpr_train, color='purple', lw=2, label=f'final (train) \nAUC = {roc_auc_train:.2f}')
        if baseline is not None:
            fpr_b, tpr_b, _ = roc_curve(y_test, baseline)
            baseline_roc = auc(fpr_b, tpr_b)
            plt.plot(fpr_b, tpr_b, color='red', lw=2, label=f'baseline (test) \nAUC = {baseline_roc:.2f}')
        plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0]); plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
        plt.title('Performance on Test Set'); plt.legend(loc='lower right')
        plt.savefig(self.out('roc.pdf'))
        plt.close()
        # plt.show()

    def plot_pr(self, y_test, proba, baseline=None, y_train=None, proba2=None):
        """
        Save Precision-Recall curves for:
        - test predictions of the evaluated model,
        - (optional) train predictions of the evaluated model,
        - (optional) baseline model predictions on test data.

        Parameters
        ----------
        y_test : pandas.Series or numpy.ndarray
            True labels for the test set.
        proba : numpy.ndarray
            Predicted probabilities for the positive class on X_test.
        baseline : numpy.ndarray or None, default None
            Baseline model probabilities on the same test labels, if available.
        y_train : pandas.Series or numpy.ndarray or None, default None
            Training labels (for plotting train PR).
        proba2 : numpy.ndarray or None, default None
            Predicted probabilities on X_train for the evaluated model.

        Returns
        -------
        None
            Saves 'pr.pdf' to figures directory.

        """
        precision, recall, _ = precision_recall_curve(y_test, proba)
        pr_auc = auc(recall, precision)
        if y_train is not None and proba2 is not None:
            precision_train, recall_train, _ = precision_recall_curve(y_train, proba2)
            pr_auc_train = auc(recall_train, precision_train)

        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        plt.figure(figsize=(10, 10))
        plt.plot(recall, precision, color='blue', lw=2, label=f'final (test) \nAUC = {pr_auc:.2f}')
        if y_train is not None and proba2 is not None:
            plt.plot(recall_train, precision_train, color='purple', lw=2, label=f'final (train) \nAUC = {pr_auc_train:.2f}')
        if baseline is not None:
            precision_b, recall_b, _ = precision_recall_curve(y_test, baseline)
            baseline_auc = auc(recall_b, precision_b)
            plt.plot(recall_b, precision_b, color='red', lw=2, label=f'baseline (test) \nAUC = {baseline_auc:.2f}')
        plt.xlim([0.0, 1.0]); plt.ylim([0.0, 1.05])
        plt.xlabel('Recall'); plt.ylabel('Precision')
        plt.title('Performance on Test Set'); plt.legend(loc='lower right')
        plt.savefig(self.out('pr.pdf'))
        plt.close()
        # plt.show()

    def plot_top_features(self, df_importances_weight, str_out):
        """
        Save bar plot for top features by a specified importance column (e.g., Gain
        or Weight). 

        Parameters
        ----------
        df_importances_weight : pandas.DataFrame
            Two-column DataFrame with importance metric and 'Feature'.
        str_out : str
            Label to include in the output file name (e.g., 'Gain_' or 'Weight_').

        Returns
        -------
        None
            Saves 'top10_features_by_{str_out}.pdf'.

        """
        cols = df_importances_weight.columns
        df_weight_sort = df_importances_weight.sort_values(cols[0], ascending=False)
        plt.rcParams['axes.labelsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.rcParams['xtick.labelsize'] = '20'
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.figure(figsize=(10, 6))
        plt.bar(x=df_weight_sort[cols[1]], height=df_weight_sort[cols[0]])
        plt.ylabel('Top 10 / ' + str(len(df_importances_weight)) + f' {cols[1]}s', rotation=90)
        plt.xlabel(f'Final Model {cols[0]}s (F-Score)')
        plt.xticks(rotation=90, ha='right')
        plt.yticks(rotation=90)
        plt.savefig(self.out(f'top10_features_by_{str_out}.pdf'))
        plt.close()
        # plt.show()

    def plot_confusion(self, y_true, proba, name, boundary):
        """
        Save  a confusion matrix plot using a specified decision threshold.

        Parameters
        ----------
        y_true : pandas.Series or numpy.ndarray
            Ground-truth binary labels.
        proba : numpy.ndarray
            Predicted probabilities for the positive class.
        name : str
            Tag for the plot filename (e.g., 'test', 'train', 'all').
        boundary : float, default 0.5
            Decision threshold used to binarize probabilities.

        Returns
        -------
        None
            Saves 'confusion_{name}_{boundary}.pdf'.

        """
        prediction = (proba >= boundary).astype(int)
        cm = confusion_matrix(y_true, prediction)
        plt.rcParams['axes.labelsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.rcParams['xtick.labelsize'] = '20'
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.figure(figsize=(10, 10))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm)
        disp.plot(cmap=plt.cm.Blues)
        plt.title(f'model boundary at {int(boundary*100)}')
        plt.savefig(self.out(f'confusion_{name}_{round(boundary,2)}.pdf'))
        plt.close()
        # plt.show()

    def plot_f1_and_pick_threshold(self, y_true, proba, xline) -> float:
        """
        Sweep decision thresholds in [0, 1] to compute F1 scores, plot the
        F1 curve, and return the threshold that maximizes F1.

        Parameters
        ----------
        y_true : pandas.Series or numpy.ndarray
            Ground-truth binary labels.
        proba : numpy.ndarray
            Predicted probabilities for the positive class.
        xline : float, default 0.5
            Vertical reference line to draw on the plot.

        Returns
        -------
        float
            The threshold value that gives maximum F1 score.

        """
        thresholds = np.linspace(0, 1, 100)
        f1_scores = []
        for t in thresholds:
            y_pred = (proba >= t).astype(int)
            f1_scores.append(f1_score(y_true, y_pred))

        # get max_f1, with highest threshold (tie-breaker)
        max_f1 = max(f1_scores)
        candidate_thresholds = thresholds[np.isclose(f1_scores, max_f1)]
        max_f1_threshold = float(candidate_thresholds.max())

        plt.rcParams['axes.labelsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.rcParams['xtick.labelsize'] = '20'
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.figure(figsize=(10, 10))
        plt.plot(thresholds, f1_scores, marker='o')
        plt.xlabel('Model Boundary')
        plt.ylabel('F1 Score')
        plt.title('F1 Score vs. Model Boundary')
        plt.axvline(x=max_f1_threshold)
        plt.axvline(x=xline, c='red')
        plt.grid(True)
        plt.savefig(self.out('f1thresh.pdf'))
        plt.close()
        # plt.show()
        self.threshold = max_f1_threshold
        return max_f1_threshold

    def plot_classification_report(self, report, threshold):
        """
        Save bar plot summarizing precision, recall, and accuracy for baseline/test/train
        (as prepared by `gbdt_classify`).

        Parameters
        ----------
        report : pandas.DataFrame
            DataFrame with rows as categories and columns as metrics, already
            aggregated/scaled to percentages.
        threshold : int, default 50
            Threshold (as a percent) to display in the plot title.

        Returns
        -------
        None
            Saves 'classif_report.pdf' to figures directory.

        """
        plt.rcParams['axes.labelsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.rcParams['xtick.labelsize'] = '20'
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.figure(figsize=(10, 6))
        report.T.plot(kind='bar', title=f'Classification Report ({threshold}%)')
        plt.ylabel('%')
        plt.xticks(rotation=90)
        plt.ylim([0, 100])
        plt.savefig(self.out('classif_report.pdf'))
        plt.close()
        # plt.show()

def get_relative_paths()->tuple[str, str, str]:
    '''
    Helper function to get relative paths for ML scripts.
    
    '''
    path = os.path.join(os.path.dirname(__file__), '..', 'data')
    model_out = os.path.join(path, '..' ,'scoring_models')
    dpath = os.path.dirname(__file__)
    return path, model_out, dpath

def clean_drug_index(df)->pd.DataFrame:
    '''
    Helper function to format DE matrix index when loading.
    
    '''
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

def get_inputs(path)->tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    '''
    Load differential expression profiles and proteome from HBD screen.
    
    '''
    LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 
                                                        'Drug_LFCxadjPval_250305a.csv'
                                                        ), delimiter = ';', decimal = ',', index_col=0, header = 0).T)
    expression_matrix = pd.read_csv(os.path.join(path,
                                                'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'
                                                ),delimiter = ',', decimal = '.', index_col=0, header = 0).T
    
    # Load metadata containing unformatted response var
    AZmeta= pd.read_csv(os.path.join(path, 
                                     'AZcompound_metadata_clustered_240611a.tsv'
                                     ),index_col=0)
    
    return LFC_matrix, expression_matrix, AZmeta

def format_response_var(AZmeta, LFC_matrix)->pd.DataFrame:
    '''
    Formats response variable and aligns to predictor index.
    
    '''
    AZmeta.index = AZmeta.index.str.replace('-','_') # Format index
    AZmeta['Gal'] = AZmeta['Gal'].str.replace('>','').astype(float) # Convert IC50 to float
    IC50s = pd.DataFrame({'Gal_IC50': AZmeta['Gal']}, dtype = float) # Extract IC50s
    IC50s = IC50s.loc[LFC_matrix.index] # Match to proteome index
    IC50snoNA = IC50s[IC50s['Gal_IC50'].isna()==False] # remove NAs
    BinaryTox = pd.DataFrame({ 'IC50' : (IC50snoNA['Gal_IC50']<10).astype(int)}) # binarize
    return BinaryTox

def extract_genes(df, pathways, matrix, boundary, outpath)->list[str]:
    """
    For each pathway enriched in split, the function collects the listed proteins, 
    and saves the final set to 'enriched_proteins.csv'.

    Parameters
    ----------
    df : pandas.DataFrame
    pathways : list of str
    matrix : pandas.DataFrame
    boundary : float
    outpath : str

    Returns
    -------
    List of genes from enrichment analysis on split. 

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
    pd.Series(selected_genes).to_csv(os.path.join(outpath, 'enriched_proteins.csv'), index=0)
    return list(selected_genes)
