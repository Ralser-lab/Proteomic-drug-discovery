# %% Establish absolute path
import sys
import os
import argparse
sys.path.append(os.path.abspath("src"))

# Import devices
from device_derivativefilter import DerivativeFilter
from device_imputer import Imputer
from device_batcher import Batcher
from device_summarizer import Summarizer

# For read into pipeline
import pandas as pd

# %% Pipeline Assembly

def main(matrix_fname, batchdata_fname, output_fname):

    dire = os.path.abspath(os.path.join(os.getcwd())) # path to repo 

    directory = os.path.join(dire, 'input') # path for input folder

    directory2 = os.path.join(dire, 'output') # path for output folder

    # Failsafe if no arguments are parsed, works on local repo / unittest folder

    if matrix_fname is None:

        matrix_fname = 'SB_PROTAC_prmatrix_plateswap_240606a.tsv'

    if batchdata_fname is None:

        batchdata_fname = '20240314_AF_50-0121_metadata_plateSwapRequested.xlsx'

    if output_fname is None:

        output_fname = 'SB_PROTAC_prmatrix'

    print(f"Matrix file: {matrix_fname}")

    print(f"Batch data file: {batchdata_fname}")

    print(f"Output file prefix: {output_fname}")

    # Load MS fragment intensities

    input_data = pd.read_csv(os.path.join(directory, matrix_fname),
                             delimiter='\t', 
                             index_col=[0,1], 
                             header=[0]).T   

    # Filter frags based on detection

    input_filter = DerivativeFilter(input_data) 

    input_sample_stats = input_filter.calculate_sample_stats(step=20)

    input_filtered = input_filter.apply_sample_filter(input_sample_stats, ramp=95)

    input_filter.save_output(input_filtered.T, os.path.join(directory2, output_fname + '_filtered_95'))

    # Impute missing frags

    input_imputer = Imputer(input_filtered, threshold=50)

    input_imputer.detection_probability_curve(boundary=0.5)

    input_imputer.precursor_missing_matrix(plot=True)

    input_imputed = input_imputer.impute_missing_matrix(knn=3)

    input_imputer.save_output(input_imputed.T, os.path.join(directory2, output_fname + '_filtered_95_imputed_50'))

    # Batch-correct on MS run date

    input_batcher = Batcher(input_imputed, 
                            path=os.path.join(directory, batchdata_fname),
                            batchID='MS.Batch',
                            logtransform=True)
    
    input_batcher.batch_correct(toPlot=True)

    input_batcher.save_output(input_batcher.input.T, os.path.join(directory2, output_fname + '_filtered_95_imputed_50_ltrfm'))
    
    input_batcher.save_output(input_batcher.output.T, os.path.join(directory2, output_fname + '_filtered_95_imputed_50_ltrfm_batched'))

    # Summarize to final proteome

    input_summarizer = Summarizer(input_batcher.output, 
                                  lfqscript=os.path.join(dire, 'maxLFQ.R'),
                                  directory=dire)
    
    input_batched_long = input_summarizer.input_long

    path = input_summarizer.save_output(input_batched_long, os.path.join(directory2, output_fname + '_filtered_95_imputed_50_ltrfm_batched_long'))

    input_summarizer.maxlfq(longform=path, convert=False) # Set to true if you need to swap from UNIPROT to SYMBOL

# %% For commandline run

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some files.')

    parser.add_argument('-m', '--matrix', help='Matrix file name', default=None)

    parser.add_argument('-b', '--batchdata', help='Batch data file name', default=None)

    parser.add_argument('-o', '--output', help='Output file name prefix', default=None)

    args = parser.parse_args()

    main(args.matrix, args.batchdata, args.output)