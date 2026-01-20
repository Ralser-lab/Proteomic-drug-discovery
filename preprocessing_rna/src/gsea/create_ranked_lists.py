# %%
import os
import pandas as pd
import numpy as np
from config_loader import load_config
from loguru import logger

# %%

def split_contrast_lists(
        data_dir: str,
        bh_path:str,
        lfc_path:str,
        input_dir:str='input',
    )->None:
    '''
    Function that processes limma outputs to create contrast lists, 
    that are rank ordered by sign(LFC) * -log10(adjusted-pvalue). 

    Contrast lists are formatted ready for input by gseapy.prerank, and written
    to disk 

    Parameters:
    ----------
    data_dir : str
        Path to directory where limma & lfc output files
    bh_path : str
        Path to limma output file containing BH adjusted p-values.
    lfc_path : str
        Path to lfc output file containing Fold Change values.
    input_dir : str
        Path to write individual contrast lists, default = 'input' in cwd().
    
    '''

    # Relative path.
    rel_path =  os.path.dirname(__file__)

    logger.info('Reading bh p-val & lfc dfs')

    # Create singed(lfc) x -log10(bh-adjusted p-value) matrix.
    bh_df = pd.read_csv(os.path.join(
                                    rel_path, # from function parameters
                                    '..',
                                    '..',
                                    data_dir, # from function parameters
                                    bh_path, # from function parameters
                                     ),
                                       index_col = 0)

    lfc_df = pd.read_csv(os.path.join(
                                    rel_path, # from function parameters
                                    '..',
                                    '..',
                                    data_dir, # from function parameters
                                    lfc_path, # from function parameters
                                     ),
                                       index_col = 0)
    
    logger.info('Combining into -log10(bh-adjusted pval) * signed(lfc) matrix')
    lfc_bh_df = np.sign(lfc_df) * -np.log10(bh_df)

    # Clean up index labels.
    lfc_bh_df.index = lfc_bh_df.index.str.replace('group', '')

    # Write sorted lists to disk. 
    for i in lfc_bh_df.index:
        logger.info(f"Writing {i} contrast to disk...")
        lfc_bh_df.loc[i].sort_values().to_csv(
            os.path.join(
                rel_path, 
                '..',
                input_dir, # from function parameters
                f'{i.replace(" ","_")}.csv' # remove spaces from filename
            ))
        logger.info(f"Writing done")

    return 

if __name__ == '__main__':
    config = load_config()

    # Send logs to a file
    logger.add(f"{config['log_dir']}/create_ranked_lists.log", 
            format="{time} {level} {message}", 
            level="INFO")

    split_contrast_lists(
        data_dir = config['data_dir'],
        lfc_path=config['lfc_path'],
        bh_path = config['bh_path'],
        input_dir = config['input_dir'])
