# DIA-MS-tools

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/docs/pipe.png)

## Description

Data-driven pipeline to build complete proteomes with reduced batch effects using precursor fragment intensities detected by data-independent-acquisition mass-spectrometry (DIA-MS). Resulting abundance follows a normal probability distribution. 

Tested on 1000, 2000, and 5000 sample drug screens so far (subjected to DIA-MS [dia-PASEF & scanning SWATH modes] and DIA-NN) - more samples (1000+) & fragments (500000+) the better.

Modules include QC, filtration, imputation, batch correction, and a wrapper for maxLFQ algorithm.

## Requirements

- Python 3.11.5 (download dependencies in local or virtual environment from 'requirements.txt') 

- R 4.3.1 (downloads required packages from CRAN at run time)

## Installation 

Clone 'DIA-MS-tools' locally and cd into root directory. Before running, make an input subfolder (named 'input') and copy the following two files:

### Precursor matrix (.tsv)

Precursor matrix containing sample IDs as columns and precursor intensities as rows. Rows should be indexed by protein and fragment ID (first and second columns). 

### Batch vector (.tsv)

Metadata file which contains batch IDs as a categorical feature (column) and sample IDs as rows. Target column must be named 'MS.Batch'.

## Execution

Create an output subfolder in directory (named 'output'), and execute main.py.

python3 main.py -m <matrix_file> -b <batch_data_file> -o <output_prefix>

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/docs/screen.png)

