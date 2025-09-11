# %% Most imports for all the devices
from datetime import datetime
import pandas as pd
import numpy as np
import math 

# Device with global functions (CV calc, plotting, saving, ...) to inherit
class Device: 

    def __init__(self):

        return

    def round_up(self, n, decimals=0):

        multiplier = 10**decimals

        return math.ceil(n * multiplier) / multiplier

    def calculate_cv(self, df):

        # Calculate CV = (standard deviation / mean) * 100

        # Replace 0 with NaN to avoid division by zero errors

        means = df.mean().replace(0, np.nan)

        stds = df.std()

        cv = (stds / means) * 100

        return cv
    
    def save_output(self, df, filename):
        
        tosave = df.copy()

        today_date = datetime.now()

        formatted_date = today_date.strftime("%y%m%d")

        path = filename + '_' + formatted_date + '.tsv'

        filename_dated  = path.split('/')[-1]

        tosave.to_csv(path, sep = '\t', index = True)

        return filename_dated
    
    def calculate_cv_ID(self, frame, name):

        means = frame.mean(axis=0)

        stdev = frame.std(axis=0)

        plot = pd.DataFrame({'Means':means,
                         'Stdev':stdev,
                         'ID':name})
        
        return plot