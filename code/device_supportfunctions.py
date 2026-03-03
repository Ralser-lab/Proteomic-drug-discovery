import os
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import pandas as pd
from sklearn.metrics import precision_recall_curve, average_precision_score

class GBDTUtils:
    """
    Utility functions that work outside of GBDT model class.
    
    """
    @staticmethod
    def cv_perform(round1, round2, workflow):
        """
        Compare performance across CV folds of two models with boxplot. 

        """
        scores1 = round1.best_cv_scores_
        mean1   = round1.best_cv_mean_
        sd1     = round1.best_cv_sd_

        scores2 = round2.best_cv_scores_
        mean2   = round2.best_cv_mean_
        sd2     = round2.best_cv_sd_

        print("Round1 fold scores:", scores1)
        print(f"Round1 mean ± SD: {mean1:.4f} ± {sd1:.4f}")
        print()
        print("Round2 fold scores:", scores2)
        print(f"Round2 mean ± SD: {mean2:.4f} ± {sd2:.4f}")

        plt.figure(figsize=(3,4))
        plt.boxplot([scores1, scores2], widths=0.5)
        # jittered fold points
        x1 = np.random.normal(1, 0.03, size=len(scores1))
        x2 = np.random.normal(2, 0.03, size=len(scores2))

        plt.scatter(x1, scores1)
        plt.scatter(x2, scores2)

        plt.ylim(0,1)
        plt.ylabel("Average Precision", fontsize=14)
        plt.xticks([1, 2], ["Round1", "Round2"], fontsize=12)
        plt.yticks(fontsize=12)

        plt.title(
            f"5-fold CV Performance\n"
            f"R1: {mean1:.3f} ± {sd1:.3f}   "
            f"R2: {mean2:.3f} ± {sd2:.3f}",
            fontsize=14
        )

        plt.tight_layout()
        plt.savefig(os.path.join(
            os.path.dirname(__file__), '..', 'figures', 
            f'{workflow}_cv_performance_plot.pdf')
        )

    @staticmethod
    def plot_pr_two_models(
        model1,
        model2,
        X_test1,
        X_test2,
        y_test,
        workflow,
        labels=("Round1", "Round2"),
        ylim=(0.0, 1.05)
    ):
        """
        Plot Precision–Recall curves for two fitted models on the same test labels.

        Parameters
        ----------
        model1, model2 : fitted estimators
            Must implement predict_proba.
        X_test1, X_test2 : array-like
            Test features for model1 and model2 (can differ if feature sets differ).
        y_test : array-like
            True binary labels (0/1).
        workflow : str
            Suffix to save figure.
        labels : tuple[str, str]
            Legend labels for model1 and model2.
        filename : str
            Output file name.
        ylim : tuple[float, float]
            y-axis limits for precision.
        """

        file_out = os.path.join(os.path.dirname(__file__), '..', 'figures')

        def _pos_proba(model, X):
            proba = model.predict_proba(X)
            if proba.ndim != 2 or proba.shape[1] < 2:
                raise ValueError("predict_proba must return shape (n_samples, 2+) for binary PR curve.")
            return proba[:, 1]  # positive-class probability

        # probabilities
        p1 = _pos_proba(model1, X_test1)
        p2 = _pos_proba(model2, X_test2)

        # PR curve points + AP
        prec1, rec1, _ = precision_recall_curve(y_test, p1)
        prec2, rec2, _ = precision_recall_curve(y_test, p2)

        ap1 = average_precision_score(y_test, p1)
        ap2 = average_precision_score(y_test, p2)

        # baseline = prevalence
        prevalence = float(np.mean(y_test))

        plt.figure(figsize=(10, 10))

        plt.plot(rec1, prec1, lw=2, label=f"{labels[0]} (AP = {ap1:.2f})")
        plt.plot(rec2, prec2, lw=2, label=f"{labels[1]} (AP = {ap2:.2f})")

        # baseline line
        plt.hlines(prevalence, 0, 1, linestyles="--", lw=2, label=f"Baseline (AP ≈ {prevalence:.2f})")

        plt.xlim(0.0, 1.0)
        plt.ylim(*ylim)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision–Recall Curve (Test Set)")
        plt.legend(loc="lower left")

        plt.tight_layout()
        plt.savefig(os.path.join(file_out, f'{workflow}_pr_2_rounds.pdf'))
        plt.close()

@dataclass
class GBDTDataLoader:
    '''
    Prepared inputs for a GBDT model.
    
    '''
    path : str
    model_out: str
    dpath: str
    AZmeta: pd.DataFrame
    X: pd.DataFrame
    y: pd.Series  
    path_genes: list[str]
    
def get_relative_paths()->tuple[str, str, str]:
    '''
    Helper function to get relative paths for ML scripts.
    
    '''
    path = os.path.join(os.path.dirname(__file__), '..', 'data')
    model_out = os.path.join(path, '..' ,'scoring_models')
    dpath = os.path.dirname(__file__)
    return path, model_out, dpath
    
def get_inputs(path)->tuple[pd.DataFrame, pd.DataFrame, list[str]]:
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
    
    # Extract features from enrichment analysis in split and write to disk 
    enrichment_profiles_in_split = pd.read_csv(os.path.join(path, 'top10_FDR_reactome.csv'))
    enriched_pathways_in_split = list(enrichment_profiles_in_split.iloc[:,0].value_counts().index)
    proteins = extract_proteins(enrichment_profiles_in_split, enriched_pathways_in_split, expression_matrix, 11.8, path)

    return LFC_matrix, AZmeta, proteins

def clean_drug_index(df)->pd.DataFrame:
    '''
    Helper function to format DE matrix index when loading.
    
    '''
    df.index = df.index.map(lambda x: '_'.join(x.split('_')[1:3]) if len(x.split('_')) > 2 else None)
    df.index = df.index.str.replace(' - Compound', '')
    return df

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

def extract_proteins(df, pathways, matrix, boundary, outpath)->list[str]:
    '''
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

    '''
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


def load_gbdt_inputs() -> GBDTDataLoader:
    '''
    Loads inputs, formats y, aligns X and y, and returns a structured dataset object.

    '''
    path, model_out, dpath = get_relative_paths()
    LFC_matrix, AZmeta, path_genes = get_inputs(path)
    y = format_response_var(AZmeta, LFC_matrix)
    X = LFC_matrix.loc[y.index]

    return GBDTDataLoader(
        path=path,
        model_out=model_out,
        dpath=dpath,
        AZmeta=AZmeta,
        X=X,
        y=y,
        path_genes=path_genes
    )

