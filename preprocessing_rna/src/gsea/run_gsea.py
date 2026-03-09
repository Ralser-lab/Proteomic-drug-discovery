# %%
import os
import pandas as pd
import gseapy as gp
from gseapy.gsea import Prerank
from config_loader import load_config
import argparse
from loguru import logger
import time
import functools

# %%

def log_time(func):
    '''
    Timing decorator.
    '''
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        logger.info(f"{func.__name__} took {elapsed:.3f} seconds")
        return result
    return wrapper

@log_time
def gsea(
        rankingfile: pd.Series,
        gene_set: str,
    ) -> Prerank:
    """
    Run GSEA (preranked).

    Parameters:
    ----------
    rankingfile : pd.Series
        A pandas Series containing the ranked genes.
    gene_set : str
        Path to the gene set file in GMT format.

    Returns:
    -------
    Prerank
        GSEA prerank results object.
    """
    return gp.prerank(
        rnk=rankingfile,
        gene_sets=gene_set,
        threads=4,  # performs doesn't scale after 4-8. 
        outdir=None, 
    )

# %%

def main(ranked_path:str): 

    config = load_config()

    gs_name = config["gene_sets"].split('/')[-1].split('.')[0] 

    contrast = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            config['input_dir'],
            ranked_path
        ),
        index_col=0,
        header=0,
    ).squeeze()  # type: pd.Series

    logger.info(f"Starting {argparser.parse_args().c} GSEA analysis...")
    # Run GSEA
    gsea_obj = gsea(
        rankingfile=contrast,
        gene_set=config["gene_sets"], # gene sets path from config
    )

    logger.info("GSEA analysis complete, writing results to disk...")
    gsea_obj.res2d.to_csv(
        os.path.join(
            os.path.dirname(__file__),
            config["output_dir"], # output directory from config
            f'{argparser.parse_args().c.split("/")[-1].replace(".csv", "")}_gsea_results_{gs_name}.csv',
        )
        )

if __name__ == "__main__":

    # Parse CLI arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-c",
        type=str,
        required=True,
        help="Path to the ranking file for GSEA.",
    )

    ranked_path = argparser.parse_args().c

    main(ranked_path)

