import os
import pandas as pd

def get_relative_paths()->tuple[str, str, str]:
    '''
    Helper function to get relative paths for ML scripts.
    
    '''
    path = os.path.join(os.path.dirname(__file__), '..', 'data')
    model_out = os.path.join(path, '..' ,'scoring_models')
    dpath = os.path.dirname(__file__)
    return path, model_out, dpath

def clean_drug_index(df)->pd.DataFrame:
    '''
    Helper function to format DE matrix index when loading.
    
    '''
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

def get_inputs(path)->tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    '''
    Load differential expression profiles and proteome from HBD screen.
    
    '''
    LFC_matrix = clean_drug_index(pd.read_csv(os.path.join(path, 
                                                        'Drug_LFCxadjPval_250305a.csv'
                                                        ), delimiter = ';', decimal = ',', index_col=0, header = 0).T)
    expression_matrix = pd.read_csv(os.path.join(path,
                                                'SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_forlimma_240611a.tsv'
                                                ),delimiter = ',', decimal = '.', index_col=0, header = 0).T
    
    # Load metadata containing unformatted response var
    AZmeta= pd.read_csv(os.path.join(path, 
                                     'AZcompound_metadata_clustered_240611a.tsv'
                                     ),index_col=0)
    
    return LFC_matrix, expression_matrix, AZmeta

def format_response_var(AZmeta, LFC_matrix)->pd.DataFrame:
    '''
    Formats response variable and aligns to predictor index.
    
    '''
    AZmeta.index = AZmeta.index.str.replace('-','_') # Format index
    AZmeta['Gal'] = AZmeta['Gal'].str.replace('>','').astype(float) # Convert IC50 to float
    IC50s = pd.DataFrame({'Gal_IC50': AZmeta['Gal']}, dtype = float) # Extract IC50s
    IC50s = IC50s.loc[LFC_matrix.index] # Match to proteome index
    IC50snoNA = IC50s[IC50s['Gal_IC50'].isna()==False] # remove NAs
    BinaryTox = pd.DataFrame({ 'IC50' : (IC50snoNA['Gal_IC50']<10).astype(int)}) # binarize
    return BinaryTox

def extract_genes(df, pathways, matrix, boundary, outpath)->list[str]:
    """
    For each pathway enriched in split, the function collects the listed proteins, 
    and saves the final set to 'enriched_proteins.csv'.

    Parameters
    ----------
    df : pandas.DataFrame
    pathways : list of str
    matrix : pandas.DataFrame
    boundary : float
    outpath : str

    Returns
    -------
    List of genes from enrichment analysis on split. 

    """
    all_genes = set()
    for pathway in pathways:
        path_genes = df[df[
            'term description'] == pathway]['matching proteins in your input (labels)'
                                            ].values[0].split(',')
        all_genes.update(set(path_genes))
    mask = matrix.mean(axis=0) > boundary
    strong_genes = matrix.columns[mask]
    selected_genes = [g for g in matrix.columns if g in strong_genes and g in all_genes]
    pd.Series(selected_genes).to_csv(os.path.join(outpath, 'enriched_proteins.csv'), index=0)
    return list(selected_genes)
