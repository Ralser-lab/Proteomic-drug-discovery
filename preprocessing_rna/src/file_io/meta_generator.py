# %%
import pandas as pd

def meta_gen(data_df:pd.DataFrame) -> pd.DataFrame:
    """
    Create a metadata table by parsing sample IDs from the index of a DataFrame 
    that contains:
        - rep        (first field)
        - cell_line  (second field, with C4_2 → C4-2)
        - drug       (third field)
        - dose       (second-to-last field, with 'Vehicle' → '0')
        - time       (last field)

    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame whose index contains sample IDs.

    Returns
    -------
    pd.DataFrame
        Metadata
    """
    return pd.DataFrame({"id": data_df.index.astype(str).str.replace('C4_2', 'C4-2')},
                index=data_df.index).assign(
        rep = lambda df: df['id'].str.split('_').str[0],
        cell_line = lambda df: df['id'].str.split('_').str[1],
        drug = lambda df: df['id'].str.split('_').str[2],
        dose = lambda df: df['id'].str.split('_').str[-2].str.replace('Vehicle','0'),
        time = lambda df: df['id'].str.split('_').str[-1],

    )
if __name__ == '__main__':
    data_df = pd.read_csv('../../data/raw_counts.csv', index_col = 0).T
    meta_gen(data_df).to_csv('../../data/metadata.csv')
