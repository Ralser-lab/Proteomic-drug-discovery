# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device import Device

# third party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import seaborn as sns
from sklearn.impute import KNNImputer
import statsmodels.api as sm

class Imputer(Device):
   
    def __init__(self, sample_filtered, threshold):

        self.data = sample_filtered

        self.boundary = None

        self.means = None

        self.means_log2 = None

        self.missing_matrix = None

        self.ramp = threshold

        self.missing = None

        self.plot_matrix = None

        self.force = None

        return
    
    def detection_probability_curve(self, p = 0.001, boundary = 0.5):

        ramp = self.ramp

        # drop introduced precursor that are missing in all samples, that were introduced after sample filtration 

        print('Remove ' + str(self.data.isna().all().sum()) + ' features with 100% NaN values / sample')
        
        self.data = self.data.dropna(axis = 1, how = 'all')

        # set 0 values in pr matrix to NaN

        print('Replace ' + str((self.data == 0).sum().sum()) + '  features with 100% zero abundance values / sample with NaN')

        self.data = self.data.replace(0, np.nan)

        # calculate precursor stats for DPC curve

        self.means = self.data.mean(axis=0)

        self.means_log2 = np.log2(self.means.astype(float))

        self.means_log2.replace([np.inf, -np.inf], np.nan, inplace=True)

        self.missing = self.data.isnull().mean(axis=0)

        missing_df = pd.DataFrame({
            'Avg Precursor Abundance' : self.means_log2.values, 
            'Missing %' : 1-(self.missing.values)})

        # Remove NA, Fit model

        cleaned_df = missing_df.dropna(subset=['Avg Precursor Abundance'])

        self.modelinput = missing_df

        X = cleaned_df[['Avg Precursor Abundance']]

        y = cleaned_df['Missing %']
             
        # Add a constant to the model for the intercept
        X = sm.add_constant(X)

        # Fit the Logistic Regression Model using GLM
        model = sm.GLM(y, X, family=sm.families.Binomial(sm.families.links.logit()))

        result = model.fit()

        # Extract the intercept and coefficient
        intercept, coef = result.params

        # Set the desired probability
        # probability / boundary = 0.4

        # Solve for the value of 'Avg Precursor Abundance' that gives y = 40%
        # log(odds) = ln(p / (1 - p)) -> Solve for 'Avg Precursor Abundance'
        # ln(0.4 / 0.6) = intercept + coef * Avg_precursor_abundance

        avg_precursor_abundance = (np.log(boundary / (1 - boundary)) - intercept) / coef

        self.boundary = self.round_up(avg_precursor_abundance, 2) 

        print(result.summary())

        mnar = cleaned_df[cleaned_df['Avg Precursor Abundance'] < self.boundary]

        mar = cleaned_df[cleaned_df['Avg Precursor Abundance'] >= self.boundary]

        omit = cleaned_df[cleaned_df['Missing %'] <= ramp/100]

        # Generating a range of values for Avg Precursor Abundance

        x_range = np.linspace(X['Avg Precursor Abundance'].min(), X['Avg Precursor Abundance'].max(), 100)

        x_range = sm.add_constant(pd.DataFrame(x_range, columns=['Avg Precursor Abundance']))

        # Predicting probabilities

        y_pred = result.predict(x_range)

        plt.figure(figsize=(12, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 
        
        plt.plot(x_range['Avg Precursor Abundance'], y_pred, marker = 'o', label = 'p-val: ' + str(p))  

        plt.axvline(self.boundary, 
            label = 'Decision boundary (' + str(boundary*100) + '%): ' + str(self.boundary), 
            color = 'red', linestyle = '--')
        
        plt.axhline(ramp/100,
            color = 'black', linestyle = '--')

        plt.scatter(mar['Avg Precursor Abundance'],
                    
            mar['Missing %'], alpha=0.1, s = 1, color = 'orange', label = 'knn(3)-impute')
          
        plt.scatter(mnar['Avg Precursor Abundance'],
             mnar['Missing %'], alpha=0.1, s = 1, color = 'yellow', label = 'min-impute')
        
        plt.scatter(omit['Avg Precursor Abundance'],
             omit['Missing %'], alpha=0.1, s = 1, color = 'grey', label = 'Below threshold (' + str(ramp) + '%): drop')
        
        # Labeling the plot

        plt.xlabel('Average Log2 Precursor Abundance')

        plt.ylabel('Average Precursor Decection')

        plt.title('Detection Probability Curve')   

        print('Decision boundary value (avg log2)= ' + str(self.boundary))

        plt.grid(True)

        plt.legend()

        plt.show()

        return 

    def precursor_missing_matrix(self, plot = False):

        log_abundance = self.means_log2.astype('float')

        mask_MNAR = log_abundance < self.boundary

        min_values = self.data.min()

        sample_MNAR = self.data.copy()

        sample_MNAR.loc['mask'] = mask_MNAR

        sample_MNAR.loc['min'] = min_values
        
        sample_MNAR_T = sample_MNAR.T

        subset_MNAR =  sample_MNAR_T[sample_MNAR_T['mask']==True].T

        subset_MAR = sample_MNAR_T[sample_MNAR_T['mask']==False].T

        precursor_names = subset_MNAR.columns

        subset_MNAR.columns = range(len(subset_MNAR.columns))

        subset_MNAR_fill = subset_MNAR.fillna('MNAR')

        subset_MNAR_fill.columns = precursor_names

        samples_MNAR_fill = pd.merge(subset_MAR, subset_MNAR_fill, left_index = True, right_index = True)

        samples_MAR = samples_MNAR_fill.drop(['mask','min'])

        samples_MAR_fill = samples_MAR.fillna('MAR')

        samples_MAR_fill = samples_MAR_fill.astype(str)

        samples_MAR_fill[~samples_MAR_fill.isin(['MAR','MNAR'])] = 'detected'

        self.missing_matrix  = samples_MAR_fill

        missing_matrix = self.missing_matrix

        original_columns = missing_matrix.columns.copy()

        missing_matrix.columns = range(len(original_columns))

        list = missing_matrix.eq('MAR').mean()

        list2 = missing_matrix.eq('MNAR').mean()

        false_positive_MAR = list[list > 1 - (self.ramp/100)].index

        false_positive_MNAR = list2[list2 > 1 - (self.ramp/100)].index

        missing_matrix[false_positive_MAR] = missing_matrix[false_positive_MAR].replace('MAR', 'todrop')

        missing_matrix[false_positive_MNAR] = missing_matrix[false_positive_MNAR].replace('MNAR', 'todrop')

        plot_matrix = missing_matrix.copy()

        columns_to_drop = plot_matrix.isin(['todrop']).any()

        plot_matrix = plot_matrix.drop(columns=plot_matrix.columns[columns_to_drop])

        self.plot_matrix = plot_matrix

        if plot == False: 
            
            return self.plot_matrix
        
        else:
        
            self.plot_missing_matrix()

        return self.plot_matrix

    def plot_missing_matrix(self):
            
        value_map = {'detected': 0, 'MAR': 1, 'MNAR': 2}

        numeric_matrix = self.plot_matrix.replace(value_map).astype(int)

        numeric_matrix.columns = range(len(numeric_matrix.columns))

        print(str(len(numeric_matrix.columns)) + ' precursors plotted in missing matrix')

        # color map

        cmap = ListedColormap(['purple', 'orange', 'yellow'])

        # Create legend handles

        detected_patch = mpatches.Patch(color='purple', label='detected')

        mar_patch = mpatches.Patch(color='orange', label='knn(3)-impute')

        mnar_patch = mpatches.Patch(color='yellow', label='min-impute')

        # Count the number of 'detected' (0s) in each column

        count_detected = (numeric_matrix == 0).sum()

        # Sort columns based on the count of 'detected'

        sorted_columns_by_detected = count_detected.sort_values(ascending=False).index

        # Sort the DataFrame based on the sorted columns

        numeric_matrix_sorted_by_detected = numeric_matrix[sorted_columns_by_detected]

        shape = numeric_matrix_sorted_by_detected.shape

        xtick_gap = int(self.round_up(shape[1]/50,-1))
                        
        ytick_gap = int(self.round_up(shape[0]/25,-1))

        plt.figure(figsize=(20, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        ax = sns.heatmap(numeric_matrix_sorted_by_detected, cmap=cmap, cbar=False)

        xticks = range(0, shape[1], xtick_gap)  

        yticks = range(0, shape[0], ytick_gap)

        ax.set_xticks(xticks)  # Set x-ticks positions

        ax.set_xticklabels(xticks, rotation=45)  # Set x-tick labels and rotate for readability

        ax.set_yticks(yticks)

        ax.set_yticklabels(yticks)  

        plt.title('Precursor Detection Matrix')

        plt.xlabel('Precursors')

        plt.ylabel('Samples')   

        plt.legend(handles=[detected_patch, mar_patch, mnar_patch], loc='lower right', fontsize = '20')

        plt.show()  

        return

    def impute_missing_matrix(self, knn = 3):

        threshold = self.ramp / 100

        log_abundance = self.means_log2.copy()

        mask_MNAR = log_abundance < self.boundary

        min_values = self.data.min()

        sample_MNAR = self.data.copy()

        sample_MNAR.loc['mask'] = mask_MNAR

        sample_MNAR.loc['min'] = min_values

        sample_MNAR.loc['miss_avg'] = self.missing

        sample_MNAR_T = sample_MNAR.T

        # Apply Sample Missing Threshold

        sample_MNAR_T = sample_MNAR_T[sample_MNAR_T['miss_avg'] < 1 - threshold]

        # Minimum Impute

        subset_MNAR = sample_MNAR_T[sample_MNAR_T['mask']==True].T

        subset_MAR = sample_MNAR_T[sample_MNAR_T['mask']==False].T

        precursor_names = subset_MNAR.columns.copy()

        subset_MNAR.columns = range(len(subset_MNAR.columns))

        subset_MNAR_imputed = subset_MNAR.fillna(subset_MNAR.loc['min'])

        subset_MNAR_imputed.columns = precursor_names

        samples_MNAR_imputed = pd.merge(subset_MAR, subset_MNAR_imputed, left_index=True, right_index=True)

        samples_min_imputed = samples_MNAR_imputed.drop(['mask','min','miss_avg'])

        # store columns and index before imputation

        columns = samples_min_imputed.columns 

        index = samples_min_imputed.index

        print(str(len(columns)) + ' precursor matrix subjected to mixed imputation')

        # KNN impute

        imputer = KNNImputer(n_neighbors=knn)
                             
        samples_knn3_imputed = imputer.fit_transform(samples_min_imputed)

        # return df with imputed values

        imputed_df = pd.DataFrame(samples_knn3_imputed, columns = columns, index = index, dtype = float)

        return imputed_df
    