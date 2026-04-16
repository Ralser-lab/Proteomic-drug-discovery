# Data

This directory holds all input data required to reproduce the main analysis pipeline (`run_all.sh`).

## Download

Download the [data](https://doi.org/10.6084/m9.figshare.28578113) from Figshare and place all files directly into this folder (`./data/`):


Alternatively, download automatically with `curl` (requires `curl` and `unzip`):
```bash
curl -L "https://figshare.com/ndownloader/articles/28578113/versions/1" -o data.zip
unzip data.zip -d data/
rm data.zip
```

## Contents

| File | Used by | Description |
|---|---|---|
| `SB_PROTAC_prmatrix_filtered_95_imputed_50_ltrfm_batched_summarized_240314.tsv` | protacs_01 | Processed PROTAC proteome matrix |
| `SB_PROTAC_metadata_240611a.tsv` | protacs_01, 04 | Sample metadata for HBD screen |
| `AZcompound_metadata_240611a.csv` | protacs_02 | AstraZeneca compound metadata |
| `c5.all.v2023.1.Hs.symbols.gmt` | protacs_01, fda_04 | GO gene sets (MSigDB C5) |
| `Cluster_LFCxPval_*uM_250305a.csv` | protacs_09, 10 | Signed -log10 adj. p-value matrices by concentration |
| `Cluster{1..15}_{1p0,10}uM_Limma_250305a.csv` | protacs_10 | Per-cluster limma DE results |
| `Drug_LFCxadjPval_250305a.csv` | protacs_09, 16, 17 | Drug-level signed DE matrix |
| `Figure3_*.xlsx` | protacs_22 | Wet-lab source data (Figure 3) |
| `Figure4_*.xlsx` | protacs_22 | Wet-lab source data (Figure 5) |
| `FigureED_Glu_gal_260119a.xlsx` | protacs_22 | Glucose/galactose viability data |
| `Figure4_ARdegdata_SafetyScores_250523a.xlsx` | protacs_21 | AR degradation IC50 + toxicity scores |

