# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device import Device

# third party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class DerivativeFilter(Device):

    def __init__(self, data):

        self.prmatrix = data

        self.data = None 

        self.sample_filtered = None

        self.peptide_step = None

        self.sample_step = None

        self.ramp = None

        return
    
    def calculate_peptide_stats(self, step = 20):

        summary = []

        self.peptide_step = step

        for x in range (0,step):

            gradient = x*(100/step)

            # X percent filter, peptide level

            select_peptide = self.prmatrix.isnull().mean(axis=0)*100

            peptide_filtered = self.prmatrix.loc[:,select_peptide <= gradient]

            sample_miss = peptide_filtered.isnull().mean(axis=1).mean()*100

            summary.append((peptide_filtered.shape[0], peptide_filtered.shape[1], sample_miss))

        filter_stats = pd.DataFrame(summary, columns = ['Samples','Precursors','Missing Sample Percentage']) 

        peptides = filter_stats['Precursors']

        return filter_stats
        
    def plot_peptide_stats(self, pasef, swath, label):

        swath_peptide_stats = swath

        pasef_peptide_stats =  pasef

        xmult = 100/self.peptide_step

        # Set plot parameters

        plt.figure(figsize=(12, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        plt.plot((pasef_peptide_stats.index)*xmult, pasef_peptide_stats['Precursors'],
                 label = label[0], marker = 'o')

        plt.plot((swath_peptide_stats.index)*xmult, swath_peptide_stats['Precursors'],
                 label = label[1], marker = 'o')

        plt.xlabel('Precursor Missing % Threshold')

        plt.ylabel('Recovered Precursors')

        plt.title('Precursor Threshold [Ramping 0 - 100%]')

        plt.legend()

        plt.grid(True)

        plt.show()

        return
    
    def calculate_sample_stats(self, step = 20):

        summary = []

        # Drop NA samples

        self.data = self.prmatrix.dropna(axis = 0, how = 'all')

        self.sample_step = step

        # Calculate Sample Missing Thresholds

        for x in range (0, step):

            gradient = x*(100/step)

            # X percent filter, sample level 

            select_sample = self.data.isnull().mean(axis =1)*100

            sample_filtered = self.data.loc[select_sample <= gradient,:]

            sample_miss = sample_filtered.isnull().mean(axis=0).mean()*100

            summary.append((sample_filtered.shape[0],sample_filtered.shape[1], sample_miss))

        filter_stats = pd.DataFrame(summary, columns = ['Samples','Precursors','Missing Sample Percentage'])

        return filter_stats
    
    def apply_sample_filter(self, filter_stats, ramp):

        xmult = 100 / self.sample_step

        self.ramp = ramp

        # Calculate Missing Threshold Ramp Derivatives

        samples = filter_stats['Missing Sample Percentage']

        threshold = filter_stats.index

        first_derivative = np.gradient(samples, threshold)

        second_derivative = np.gradient(first_derivative, threshold)
        
        self.stat = filter_stats

        # Set plot parameters

        plt.figure(figsize=(12, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 
    
        # Plot Precursors vs Missing Sample Percentage

        recovered = self.round_up(filter_stats.iloc[int(ramp/5)]['Samples'])

        plt.plot(threshold*xmult, filter_stats['Samples'], marker='o', label = str(recovered) + ' samples')

        plt.xlabel('Sample Missing % Threshold')

        plt.ylabel('Recovered Samples')

        plt.axhline(recovered, linestyle = '--', color = 'r', label = str(ramp) + '% ramp:\n')

        plt.axvline(ramp, color = 'r', linestyle = '--')

        plt.title('Sample Threshold [Ramping 0 - 100%]')

        plt.legend()

        plt.grid(True)

        plt.show()

        # Set plot parameters

        plt.figure(figsize=(12, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        # Plot the first derivative with a horizontal line at the inflection value

        plt.plot(threshold*xmult, first_derivative, marker='o')

        plt.xlabel('Sample Missing % Threshold')

        plt.ylabel('Rate of Change')

        plt.title('First Derivative')

        intercept = int(ramp/5)

        first_int = self.round_up(first_derivative[intercept], 2)

        second_int = self.round_up(second_derivative[intercept], 2)

        plt.axhline(first_int, color = 'r', linestyle = '--', 
            label = 'dy/dx ≈ ' + str(first_int) + ' at ' + str(ramp) + '%')
        
        plt.axvline(ramp, color = 'r', linestyle = '--')

        plt.legend()

        plt.grid(True)

        plt.show()

        # Set plot parameters

        plt.figure(figsize=(12, 10))

        plt.rcParams['axes.labelsize'] = '30'   

        plt.rcParams['axes.titlesize'] = '30' 

        plt.rcParams['legend.fontsize'] = '28'  

        plt.rcParams['xtick.labelsize'] = '26'  

        plt.rcParams['ytick.labelsize'] = '26' 

        # Plot the second derivative with a horizontal line at the inflection value

        plt.plot(threshold*xmult, second_derivative, marker='o')

        plt.xlabel('Sample Missing % Threshold')

        plt.ylabel('Acceleration')

        plt.title('Second Derivative')

        plt.axhline(second_int, color = 'r', linestyle = '--', 
            label = 'ddy/dxx ≈ ' + str(self.round_up(second_int,2)) + ' at ' + str(ramp) + '%')
        
        plt.axvline(ramp, color = 'r', linestyle = '--')

        plt.legend()

        plt.grid(True)

        plt.show()
        
        select_sample = self.data.isnull().mean(axis=1)*100

        # Filter out based on diagnosed peptdide missigness

        self.sample_filtered = self.data.loc[select_sample <= ramp,:]

        return self.sample_filtered