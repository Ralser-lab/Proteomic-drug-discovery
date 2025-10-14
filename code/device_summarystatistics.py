#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: device_summarystatistics.py
Description:
    Statistical utilities for exploratory analysis, hypothesis testing, and 
    regression-based IC50 estimation. Includes implementations of summary 
    statistics, z-tests, t-tests, one-sample tests, geometric mean, and 
    coefficient of variation.

Author: Shaon Basu
Date: 2025-09-29

Function architecture:
--------
- std_calc : Compute standard deviation using E(X^2) - E(X)^2.
- summary_stats : Report mean, standard deviation, and sample size for two series.
- geometric_mean : Calculate the geometric mean of a series.
- z_test : Perform a two-sample z-test.
- t_test : Perform Welch's two-sample t-test.
- linear_IC50 : Estimate IC50 values by linear regression on each DataFrame column.
- one_test : Perform a one-sample t-test against a population mean.
- calculate_cv : Compute coefficient of variation (CV) for a DataFrame by index.

Outputs
-------
- Console logs of test statistics and p-values
- pandas.DataFrame or pandas.Series results (depending on function)

Requirements
------------
Python >= 3.8  
Dependencies: pandas, numpy, scipy, statsmodels

"""
# %%
import os
import sys
sys.path.append(os.path.dirname(__file__))

import scipy.stats as stats # for cdf
import statsmodels.api as sm # for linear fits
import pandas as pd 
import numpy as np

def std_calc(df):
    """
    Calculate var using E(X^2) - E(X)^2 on a series.

    """
    var = (df**2).mean(skipna = True) - (df.mean(skipna = True))**2
    std = np.sqrt(var)
    return std

def summary_stats(df1, df2):
    """
    Calculate parameters (n, std, mew) for 2 sample vars.

    """
    mean1, mean2 = df1.mean(skipna = True), df2.mean(skipna = True) # calc E(X)
    std1, std2 = std_calc(df1), std_calc(df2) # calc E(X - E(X))
    n1, n2 = df1.count(), df2.count() # calc N
    print(mean1, std1, n1, mean2, std2, n2) 
    return mean1, std1, n1, mean2, std2, n2

def geometric_mean(ser1):
    """
    Calculate geometric mean of a series.

    """
    product = ser1.prod()
    return np.power(product, 1/ser1.shape[0])

def z_test(df1, df2):
    """
    Two sample z-test.

    """
    # calculate summary stats
    mean1, std1, n1, mean2, std2, n2 = summary_stats(df1, df2)
    # std for z-test
    pooled_std = np.sqrt((std1**2/n1) + (std2**2/n2)) 
    # z - score
    z_score = (mean1 - mean2) / pooled_std 
    # 2 x norm cdf(z-score) -> 2-sample z-test
    p_value = 2 * (1 - stats.norm.cdf(np.abs(z_score)))
    return z_score, p_value # return params

def t_test(df1, df2):
    """
    Two sample t-test.

    """
    # print out summary stats
    print(f"Mean and Std of df1: Mean = {df1.mean():.4f}, Std = {df1.std():.4f}")
    print(f"Mean and Std of df2: Mean = {df2.mean():.4f}, Std = {df2.std():.4f}")
    # perform independent t-test (welch by default)
    t_stat, p_value = stats.ttest_ind(df1.dropna(), df2.dropna(), equal_var =False)
    print(f'p-value: {p_value:.4f}')
    print(f't-stat: {t_stat:.4f}')
    return t_stat, p_value # return params

def linear_IC50(df):
    """
    Colwise linear fit to extract IC50, *assumes response
    is normally distributed*.

    """
    # empty vec for storing IC50s
    ic50_values = pd.Series()
    for col in df.columns:
        # drop NaNs
        valid_data = df[col].dropna()
        # dependent
        x = valid_data.index
        # response
        Y = valid_data.values
        # add intercept
        X = sm.add_constant(x)
        # regress
        model = sm.OLS(Y, X).fit()
        # extract params
        intercept, slope = model.params
        # calculate ic50
        ic50 = (0.5 - intercept)/slope
        # append series
        ic50_values[col] = ic50
    return ic50_values

def one_test(df1, popmean = 100):
    """
    One sample t-test.

    """
    # print out summary stats
    print(f"Mean and Std of df1: Mean = {df1.mean():.4f}, Std = {df1.std():.4f}")
    # perform one sample t-test using popmean
    t_stat, p_value = stats.ttest_1samp(df1.dropna(), popmean)
    print(f'p-value: {p_value:.4f}')
    print(f't-stat: {t_stat:.4f}')
    return t_stat, p_value

def calculate_cv(df, name):
    """
    Calculate coefficient of variation of a dataframe by index.
    
    """
    # CV = (standard deviation / mean) * 100 
    # Replace 0 with NaN to avoid division by zero errors
    means = df.mean().replace(0, np.nan)
    stds = df.std()
    df = pd.DataFrame({'ID': name,
                        'Means' : means,
                        'Stdev': stds})
    return df

