# %%
import rdata
import pandas as pd
import argparse
from loguru import logger 
from typing import Any

def load_rdata(path : str) -> Any:
    """
    Load an R .rds file into a Python object using the `rdata` package.

    Parameters
    ----------
    path : str
        Path to the .rds file to load.

    Returns
    -------
    Any
        A Python representation of the R object. For Bioconductor objects,
        this will typically be a `SimpleNamespace` containing nested structures.
    """
    parsed = rdata.parser.parse_file(path)
    converted = rdata.conversion.convert(parsed)
    return converted

def extract_df_rdata(r_data : dict, frametype: str = 'aligned') -> pd.DataFrame : 
    """
    Extract an assay matrix (e.g., counts, TPM, VST) from a converted 
    R SummarizedExperiment-like object and return it as a pandas DataFrame.

    Parameters
    ----------
    r_data : object
        The converted R object returned by `load_rdata`. Must contain 
        `r_data.assays.data.listData`, where assays are stored.
    frametype : str, optional
        Name of the assay to extract. 
        Common options include:
            - 'aligned'       : raw aligned read counts (integers)
            - 'counts'        : tximport estimated counts (floats)
            - 'tpm'           : transcripts per million
            - 'normalized'    : length-normalized counts
            - 'vst'           : variance-stabilized data
            - 'fpkm'          : fragments per kilobase per million
        Default is 'aligned'.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where rows are genes (dim_0) and columns are samples (dim_1).
    """
    assays_dict = r_data.assays.data.listData
    counts_da = assays_dict[frametype]
    return pd.DataFrame(
        counts_da.values,
        index=counts_da.coords["dim_0"].values,
        columns=counts_da.coords["dim_1"].values,
        )

def extract_df_tximport(tximport : dict, frametype: str = 'counts') -> pd.DataFrame : 
    """
    Convert a tximport-style assay matrix (from rdata or Python conversion)
    into a pandas DataFrame.

    Parameters
    ----------
    tximport : dict-like
        A converted R tximport list, where matrices have `.dim_0` (rownames)
        and `.dim_1` (colnames) attributes.
    frametype : str, optional
        Name of the matrix to extract, e.g.:
            - 'counts'
            - 'abundance'
            - 'length'
        Default is 'counts'.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where rows are samples and columns are genes.
        (Note: The returned matrix is transposed compared to the R version.)
    """
    return pd.DataFrame(tximport[frametype],
                            index=tximport[frametype].dim_0,
                            columns=tximport[frametype].dim_1).T
# %% Load datasets

def etl_pipeline(rds_path:str='bcb.rds')->object:

    logger.info("Starting ETL pipeline")

    logger.info("Extracting RDS...")
    r_data = load_rdata(f'../../data/{rds_path}')

    logger.info("Transforming count data to DataFrame...")
    count_df = extract_df_rdata(r_data, 'aligned')

    logger.info("Loading csv to /data...")
    count_df.to_csv('../../data/raw_counts.csv')
    logger.success("ETL pipeline completed successfully")

    return count_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", type=str, required=True, help="Path to .rds")
    args = parser.parse_args()
    etl_pipeline(args.r)