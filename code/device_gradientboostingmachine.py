# %% filepath and dependencies
###################################################################################
###################################################################################
###################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.metrics import roc_curve, auc, ConfusionMatrixDisplay, f1_score, precision_recall_curve
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
import matplotlib.pyplot as plt
import shap
import seaborn as sns
import numpy as np
import os

class GBDT: 

    def __init__(self, dpath):
        
        self.dpath = dpath    
        self.dpath = os.path.join(dpath, '..', 'data')
        self.dout = os.path.join(dpath, '..', 'figures')

    def get_model_features(self, plot = True, n = 20):

        model = self.best_model
        
        importances_weight = model.get_booster().get_score(importance_type='weight')
        importances_gain = model.get_booster().get_score(importance_type='gain')
        importances_cover = model.get_booster().get_score(importance_type='cover')

        self.df_importances = pd.DataFrame({'Weight': importances_weight, 'Gain': importances_gain, 'Cover': importances_cover}).reset_index().rename(columns={'index': 'Feature'})

        self.df_weight_sort = self.df_importances.sort_values(['Weight','Gain','Cover'], ascending = [False,False,False])[0:20]

        if plot == True:
            self.plot_top10_features(self.df_weight_sort[['Gain','Feature']], 'Gain_')
            self.plot_top10_features(self.df_weight_sort[['Weight','Feature']], 'Weight_')


        self.top_n = self.df_importances.sort_values(['Weight','Gain','Cover'],ascending = [False,False,False])[0:n]['Feature']

        return  self.top_n

    def plot_roc(self, y_test, proba, baseline = None, enriched = None, y_train = None, proba2 = None):
        # Compute ROC curve and ROC area
        fpr, tpr, thresholds = roc_curve(y_test, proba)
        roc_auc = auc(fpr, tpr)
        # Comput ROC curve and ROC area for training set, i fpassed
        if y_train is not None: 
            fpr_train, tpr_train, thresholds_train = roc_curve(y_train, proba2)
            roc_auc_train = auc(fpr_train, tpr_train)
        # Set global font size parameters using plt.rcParams
        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        # Plot ROC curve
        plt.figure(figsize=(10, 10))
        plt.plot(fpr, tpr, color='blue', lw=2, label=f'final (test) \nAUC = {roc_auc:.2f}')
        # Plot ROC curve for training set if passed
        if proba2 is not None: 
            plt.plot(fpr_train, tpr_train, color='purple', lw=2, label=f'final (train) \nAUC = {roc_auc_train:.2f}')
        # Plot ROC curve for enriched on test set if passed
        if enriched is not None:
            fpr_3, tpr_3, thresholds_3 = roc_curve(y_test, enriched)
            enriched_roc = auc(fpr_3, tpr_3)
            plt.plot(fpr_3, tpr_3, color='green',lw=2, label=f'enriched (test) \nAUC = {enriched_roc:.2f}')
            version = '3'
        # Plot ROC curve for baseline on test set if passed
        if baseline is not None:
            fpr_2, tpr_2, thresholds_2 = roc_curve(y_test, baseline)
            baseline_roc = auc(fpr_2, tpr_2)
            plt.plot(fpr_2, tpr_2, color='red',lw=2, label=f'baseline (test) \nAUC = {baseline_roc:.2f}')
            version = '2'
        else:
            version = '1'
        plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Performance on Test Set')
        plt.legend(loc='lower right')
        plt.savefig(os.path.join(self.dout, str('xgb_roc_v' + version + '.pdf')))
        plt.show()

    def plot_precisionrecall(self, y_test, proba, baseline = None, enriched = None, y_train = None, proba2 = None):
        # Compute ROC curve and ROC area
        precision, recall, thresholds = precision_recall_curve(y_test, proba)
        pr_auc = auc(recall, precision)
        # Comput ROC curve and ROC area for training set, i fpassed
        if y_train is not None: 
            precision_train, recall_train, thresholds_train = precision_recall_curve(y_train, proba2)
            pr_auc_train = auc(recall_train, precision_train)
        # Set global font size parameters using plt.rcParams
        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        # Plot ROC curve
        plt.figure(figsize=(10, 10))
        plt.plot(recall, precision, color='blue', lw=2, label=f'final (test) \nAUC = {pr_auc:.2f}')
        # Plot ROC curve for training set if passed
        if proba2 is not None: 
            plt.plot(recall_train, precision_train, color='purple', lw=2, label=f'final (train) \nAUC = {pr_auc_train:.2f}')
        # Plot ROC curve for enriched on test set if passed
        if enriched is not None:
            precision_3, recall_3, thresholds_3 = precision_recall_curve(y_test, enriched)
            enriched_auc = auc(recall_3, precision_3)
            plt.plot(recall_3, precision_3, color='green',lw=2, label=f'enriched (test) \nAUC = {enriched_auc:.2f}')
            version = '3'
        # Plot ROC curve for baseline on test set if passed
        if baseline is not None:
            precision_2, recall_2, thresholds_2 = precision_recall_curve(y_test, baseline)
            baseline_auc = auc(recall_2, precision_2)
            plt.plot(recall_2, precision_2, color='red',lw=2, label=f'baseline (test) \nAUC = {baseline_auc:.2f}')
            version = '2'
        else:
            version = '1'
        plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Performance on Test Set')
        plt.legend(loc='lower right')
        plt.savefig(os.path.join(self.dout, str('xgb_pr_v' + version + '.pdf')))
        plt.show()

    def plot_top10_features(self, df_importances_weight, str_out):
        id = df_importances_weight.columns
        df_weight_sort = df_importances_weight.sort_values(id[0], ascending = False)
        plt.rcParams['axes.labelsize'] = '20'   
        plt.rcParams['axes.titlesize'] = '20' 
        plt.rcParams['xtick.labelsize'] = '20'  
        plt.rcParams['ytick.labelsize'] = '20' 
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'

        plt.figure(figsize=(10, 6))
        plt.bar(x = df_weight_sort[id[1]], height =  df_weight_sort[id[0]]) 
        plt.ylabel('Top 10 / ' + str(len(df_importances_weight)) + f' {id[1]}s', rotation = 90)
        plt.xlabel(f'Final Model {id[0]}s (F-Score)')
        plt.xticks(rotation = 90, ha='right')
        plt.yticks(rotation = 90)
        #plt.ylim([0,13])
        plt.savefig(os.path.join(self.dout, str_out + 'xgb_top10_features.pdf'))
        plt.show()

    def plot_confusion(self, y_test, proba, name, boundary = 0.5):
        plt.rcParams['axes.labelsize'] = '20'   
        plt.rcParams['axes.titlesize'] = '20' 
        plt.rcParams['xtick.labelsize'] = '20'  
        plt.rcParams['ytick.labelsize'] = '20' 
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        prediction = (proba >= boundary).astype(int)
        cm = confusion_matrix(y_test, prediction)
        # Plot confusion matrix
        plt.figure(figsize = (10,10))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm)
        disp.plot(cmap=plt.cm.Blues)
        plt.title(f'model boundary at {boundary*100}')
        plt.savefig(os.path.join(self.dout, str('xgb_confusion_' + name + '_' + str(boundary) + '.pdf')))
        plt.show()

    def plot_f1scores(self, y_test, proba, xline = 0.5):

        thresholds = np.linspace(0, 1, 100)
        f1_scores = []

        for threshold in thresholds:
            y_pred = (proba >= threshold).astype(int)
            f1 = f1_score(y_test, y_pred)
            f1_scores.append(f1)

        # Find the index of the maximum F1 score
        max_f1_index = np.argmax(f1_scores)
        # Find the threshold that gives the maximum F1 score
        max_f1_threshold = thresholds[max_f1_index]
        plt.rcParams['axes.labelsize'] = '20'   
        plt.rcParams['axes.titlesize'] = '20' 
        plt.rcParams['xtick.labelsize'] = '20'  
        plt.rcParams['ytick.labelsize'] = '20' 
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.figure(figsize=(10, 10))
        plt.plot(thresholds, f1_scores, marker='o')
        plt.xlabel('Model Boundary')
        plt.ylabel('F1 Score')
        plt.title('F1 Score vs. Model Boundary')
        plt.axvline(x = max_f1_threshold)
        plt.axvline(x = xline, c = 'red')
        plt.grid(True)
        plt.savefig(os.path.join(self.dout, str('xgb_f1thresh.pdf')))
        plt.show()
        return max_f1_threshold

    def classification_report2(self, report, threshold = 50):
        plt.rcParams['axes.labelsize'] = '20'   
        plt.rcParams['axes.titlesize'] = '20' 
        plt.rcParams['xtick.labelsize'] = '20'  
        plt.rcParams['ytick.labelsize'] = '20' 
        plt.rcParams['ytick.labelsize'] = '20'
        plt.rcParams['legend.fontsize'] = '20'
        plt.rcParams['axes.titlesize'] = '20'
        plt.figure(figsize = (10,6))
        report.T.plot(kind = 'bar', title = f'Classification Report ({threshold}%)')
        plt.ylabel('%')
        plt.xticks(rotation = 90)
        plt.ylim([0,100])
        plt.savefig(os.path.join(self.dout, str('classif_report.pdf')))
        plt.show()

    def gbdt_baseline(self, Xb_train, Xb_test, yb_train, yb_test):

        self.Xb_train, self.Xb_test, self.yb_train, self.yb_test = Xb_train, Xb_test, yb_train, yb_test

        baseline_model = xgb.XGBClassifier(seed = 42)

        baseline_model.fit(self.Xb_train, self.yb_train)
        self.baseline_predictions = baseline_model.predict(self.Xb_test)
        self.baseline_proba = baseline_model.predict_proba(self.Xb_test)[:, 1]
        self.baseline_model = baseline_model

    def gbdt_gridcv(self, params, X_train, y_train):

        self.X_train = X_train
        self.y_train = y_train

        # Define k-fold cross-validation
        kf = KFold(n_splits=5, shuffle=True, random_state=42)

        # Perform grid search with k-fold cross-validation on the training data
        grid_search = GridSearchCV(estimator=xgb.XGBClassifier(), param_grid=params, cv=kf, scoring='accuracy', n_jobs=-1)
        grid_search.fit(X_train, y_train)

        # Get the best parameters
        best_params = grid_search.best_params_

        print(best_params)

        # Train the best model on the entire training dataset
        best_model = xgb.XGBClassifier(**best_params, seed = 42)
        best_model.fit(X_train, y_train)
        self.best_model = best_model
        return best_model
    
    def gbdt_evaluate(self, X_test, y_test, model):

        self.X_test = X_test
        self.y_test = y_test

        self.model = model

        self.proba_test = self.model.predict_proba(X_test)[:, 1]

        self.proba_train = self.model.predict_proba(self.X_train)[:,1]

        self.plot_roc(y_test = y_test, proba = self.proba_test, baseline = self.baseline_proba, 
                enriched=None, y_train = self.y_train, proba2 = self.proba_train)

        self.plot_precisionrecall(y_test = y_test, proba = self.proba_test, baseline = self.baseline_proba, 
                enriched=None, y_train = self.y_train, proba2 = self.proba_train)

    def gbdt_classify(self, X, y, model = 'best', threshold = 'auto'):

        self.X = X
        self.y = y

        if model == 'best':
            model = self.best_model
            self.proba = model.predict_proba(self.X)[:, 1]
            if threshold == 'auto':
                self.threshold = self.plot_f1scores(self.y, self.proba, xline = 0.50)
            else:
                self.plot_f1scores(self.y, self.proba, xline = threshold)
                self.threshold = threshold
        else:
            model = model
            self.proba = model.predict_proba(self.X)[:, 1]
            if threshold == 'auto':
                self.threshold = self.plot_f1scores(self.y, self.proba, xline = 0.50)
            else:
                self.plot_f1scores(self.y, self.proba, xline = threshold)
                self.threshold = threshold

        # Get predictions from best model

        self.prediction = self.proba >= self.threshold
        
        self.predict_test = (self.proba_test >= self.threshold)

        self.predict_train = (self.proba_train >= self.threshold)

        self.predict_baseline = (self.baseline_proba >= self.threshold)

        self.plot_confusion(self.y_test, self.proba_test, name = 'test', boundary = self.threshold)

        self.plot_confusion(self.y_train, self.proba_train, name = 'train', boundary = self.threshold)

        self.plot_confusion(self.y, self.proba, name = 'all', boundary = self.threshold)

        print(classification_report(self.y_test, self.predict_test))

        def format_report(df, idx):
            df['accuracy'] = df.loc['accuracy'].mean()
            df = df.iloc[:,[0,1,4]].iloc[[1]]
            df.rename(index = {'1':idx}, inplace = True)
            df = df * 100
            return df

        final_report = format_report(pd.DataFrame(classification_report(self.y_test, self.predict_test, output_dict=True)).T, 'final test')

        baseline_report = format_report(pd.DataFrame(classification_report(self.yb_test, self.predict_baseline, output_dict=True)).T, 'baseline test')

        train_report = format_report(pd.DataFrame(classification_report(self.y_train, self.predict_train, output_dict=True)).T, 'final train')

        self.report = pd.concat([baseline_report, final_report, train_report])

        self.classification_report2(self.report, round(self.threshold*100))

    def gbdt_SHAP(self, top_interactors = ['NDUFA5', 'PRKAR2B','CYC1', 'PDPR',  'MTIF2', 'PAFAH1B1','NDUFA4']):
        # SHAP values of best model

        self.top_interactors = top_interactors

        self.explainer = shap.TreeExplainer(self.best_model)
        self.shap_values = self.explainer.shap_values(self.X_test)

        # Test set, SHAP summary plot

        shap.summary_plot(self.shap_values, self.X_test, show=False, plot_size = [5, 10]) 
        plt.savefig(os.path.join(self.dout,'shap_explainer.pdf'))
        plt.close()

        # Get mean absolute SHAP values for each feature
        self.shap_abs_mean = np.abs(self.shap_values).mean(axis=0)
        
        self.shap_feature_importance = pd.DataFrame({
            'feature': self.X_test.columns,
            'importance': self.shap_abs_mean
        }).sort_values(by='importance', ascending=False)[0:20]

        #  SHAP weights for individual drugs in test set

        self.df_shap_values = pd.DataFrame(self.shap_values, index = self.X_test.index, columns = self.X_test.columns)

        # SHAP interaction plot

        interaction_values = self.explainer.shap_interaction_values(self.X_train)

        shap.summary_plot(interaction_values, self.X_train, show = True)
        plt.close()

        variances = np.var(interaction_values, axis=0)  # This calculates the variance across the samples for each feature interaction

        # Calculate variance for each interaction
        variances = np.std(interaction_values, axis=0)

        # Create DataFrame for variances
        variance_df = pd.DataFrame(variances, index=self.X_train.columns, columns=self.X_test.columns)

        # Filter the DataFrame to keep only the genes of interest (6 genes only for the 6x6 matrix, dropping "MRPS34")
        filtered_variance_df = variance_df.loc[top_interactors, top_interactors]

        # Plotting
        plt.figure(figsize=(10, 8))
        plt.rcParams['axes.labelsize'] = 35
        plt.rcParams['axes.titlesize'] = 35
        plt.rcParams['xtick.labelsize'] = 26
        plt.rcParams['ytick.labelsize'] = 26
        plt.rcParams['legend.fontsize'] = 28
        sns.clustermap(filtered_variance_df, annot=False, cmap='magma_r', linewidths=.5)
        plt.title('SHAP Interactions in Test set')
        plt.savefig(os.path.join(self.dout,'shap_interactions.pdf'))
        plt.show()
    
    def gbdt_calibrate(self):

        # Calibrate the model
        calibrated_model = CalibratedClassifierCV(self.best_model, method='isotonic', cv='prefit')  # method can be 'sigmoid' for Platt Scaling
        calibrated_model.fit(self.X_test, self.y_test)

        #Calibrate baselinemodel 
        calibrated_baseline = CalibratedClassifierCV(self.baseline_model, method = 'isotonic', cv='prefit')
        calibrated_baseline.fit(self.Xb_test, self.yb_test)
        self.baseline_predictions = calibrated_baseline.predict(self.Xb_test)
        self.baseline_proba = calibrated_baseline.predict_proba(self.Xb_test)[:, 1]

        # Predict probabilities
        prob_pos = calibrated_model.predict_proba(self.X_test)[:, 1]

        # Compute start curve
        fraction_before, mean_predicted_before = calibration_curve(self.y_test, self.proba_test, n_bins = 10)

        # Compute calibration curve
        fraction_of_positives, mean_predicted_value = calibration_curve(self.y_test, prob_pos, n_bins=10)

        # Plot calibration curve
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
        plt.savefig(os.path.join(self.dout,'calibration_curve.pdf'))
        plt.show()

        self.calibrated_model = calibrated_model

    def export_all_predictions(self, AZmeta, model):
        
        #  Look at Predictions in all data, using probabilities from calibrated model

        self.proba = model.predict_proba(self.X)[:, 1]

        self.prediction = self.proba >= self.threshold

        self.outcome = pd.DataFrame({'Actual': self.y,
                                'Predicted': self.prediction,
                                'Probability': self.proba})

        subset = self.outcome.copy()

        subset['Actual'].replace({1 : 'Toxic', 0 : 'Non-Toxic'}, inplace = True)

        subset['Drug'] = AZmeta.loc[subset.index]['Drug ID']

        subset['Cluster'] = AZmeta.loc[subset.index]['Dend']

        subset.sort_values(['Cluster'])

        # Prediction base value export for plotting in R

        clusters = subset.sort_values(['Drug','Actual'], ascending = [True, False])

        #clusters = clusters.loc[targets_kevin]

        cluster_shap = self.X.loc[clusters.index][self.shap_feature_importance['feature'].values]

        dotplotforR = pd.merge(clusters, cluster_shap, left_index=True, right_index=True)

        self.dotplotforR = dotplotforR

        dotplotforR.to_csv(os.path.join(self.dpath, 'Rplot_Figure4.csv'))

    def export_subset_predictions(self, AZmeta, targets):
        #  Look at Predictions in Test Set

        subset = self.outcome.copy()

        subset['Actual'].replace({1 : 'Toxic', 0 : 'Non-Toxic'}, inplace = True)

        subset['Drug'] = AZmeta.loc[subset.index]['Drug ID']

        subset['Cluster'] = AZmeta.loc[subset.index]['Dend']

        subset.sort_values(['Cluster'])

        # Prediction base value export for plotting in R

        clusters = subset.sort_values('Drug')

        clusters = clusters.loc[targets]

        self.get_model_features(plot = False)

        gene = list(self.df_importances.sort_values('Weight').iloc[-20:]['Feature'])

        matrix = self.X.copy()

        matrix_subset = matrix[gene]

        def plot_diffexp(column):
            sorted = column.sort_values()
            plt.rcParams['axes.labelsize'] = 35
            plt.rcParams['axes.titlesize'] = 35
            plt.rcParams['xtick.labelsize'] = 20
            plt.rcParams['ytick.labelsize'] = 26
            plt.rcParams['legend.fontsize'] = 28
            plt.figure(figsize=[10,10])
            plt.title(f'final model weight (F-1): {column.name}\n (AR PROTACs, Thalidomide 5N)')
            plt.ylabel(f'test set targets vs DMSO')
            sorted.plot(kind = 'bar')
            plt.savefig(os.path.join(self.dout, f'{column.name}_LFCxPval.pdf'))
            plt.show()

        matrix_subset.loc[targets].apply(plot_diffexp, axis = 0)

        cluster_shap = self.X.loc[clusters.index][self.shap_feature_importance['feature'].values]

        dotplotforR = pd.merge(clusters, cluster_shap, left_index=True, right_index=True)

        self.dotplotforR_subset = dotplotforR

        dotplotforR.to_csv(os.path.join(self.dpath, 'Rplot_Figure5.csv'))


# %%
