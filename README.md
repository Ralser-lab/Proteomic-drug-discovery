# Proteomic-drug-discovery

Repository containing scripts to regenerate all figures, training and analytics in 
**'Proteome-guided discovery accurately maps and mitigates toxicity mechanisms of therapeutic androgen receptor degraders'.**

## Installation 

1. Clone this repository and navigate to root directory (copy-paste in CLI):
```bash
   git clone 'https://github.com/BasuShaon/Proteomic-drug-discovery' 
   cd /Proteomic-drug-discovery
```
2. Download [data](https://figshare.com/s/6d164fd50adfdb9a68d7) and place it into `/data`.
3. Install [docker](https://www.docker.com/get-started).
4. Build the virtual environment using the pinned [Dockerfile](./docker/Dockerfile) (copy-paste in CLI):
```bash
   docker build -t prot-env -f docker/Dockerfile . 
```

## Execution

To reproduce manuscript findings, execute `run_all.sh` (copy-paste in CLI):
```bash 
   docker run --rm -v "$PWD":/image prot-env bash run_all.sh
```

## ML-modification 

To adjust hyperparameters for toxicity scoring ML workflow: 
1. Edit [`/configs/HYPER.json`](./configs/HYPER.json).
2. Then run `deepsearch.smk` in docker container (copy-paste in CLI):
```bash
   docker run --rm -v "$PWD":/image prot-env snakemake -s workflows/deepsearch.smk -j 12
```

## Environment

**Hardware** (MacBook Pro, M2 MAX 12-core CPU, 32 GB RAM, macOS Ventura 13.3) 

The following [software environment](./docker/Dockerfile) was used for development:

**Python 3.11.5** (gseapy==1.0.6, joblib==1.3.2, matplotlib==3.8.1, numpy==1.25.2, openpyxl==3.1.2, pandas==2.1.0, scikit-learn==1.3.0, scipy==1.11.2, seaborn==0.13.2, shap==0.46.0, snakemake==9.12.0, statsmodels==0.14.0, xgboost==2.0.3) 

**R 4.3.1** (ggplot2==3.5.2, dplyr==1.1.4, tidyr==1.3.1, pheatmap==1.0.13, cowplot==1.2.0, RColorBrewer==1.1-3, ggnewscale==0.5.2, ape==5.8-1, factoextra==1.0.7, ggfortify==0.4.18) 

**Bioconductor 3.18** (limma==3.58.1, EnhancedVolcano==1.20.0, ComplexHeatmap==2.18.0, ggtree==3.10.1, ggtreeExtra==1.12.0) 

## Mapping index

Generated outputs (models, figures, logfiles) save into `/scoring_models`, `/figures`, and `/logs` respectively, and they map as follows:

---

### FDA analysis (steps 1–4)

**Step 1 — `fda_01_globalvarianceanalysis.py`**
- **Outputs**
  - `figures/fda_01_dispersion_kde.pdf`
  - `figures/fda_01_global_PCA.pdf`
  - `figures/fda_01_global_PCA_scree.pdf`
  - `figures/fda_01_heatmap.pdf`
  - `figures/fda_01_heatmap.png`
  - `logs/fda_01_globalvarianceanalysis.log`

**Step 2 — `fda_02_limma_drug.R`**
- **Outputs**
  - `figures/fda_02_volcanoes_top3.png`
  - `logs/fda_02_limma_drug.log`

**Step 3 — `fda_03_de_pca.py`**
- **Outputs**
  - `figures/fda_03_2D_PCA_tmatrix.pdf`
  - `figures/fda_03_3D_PCA_tmatrix.pdf`
  - `figures/fda_03_PDF_all.pdf`
  - `figures/fda_03_PDF_stats.pdf`
  - `figures/fda_03_PMF_conditioned.pdf`
  - `figures/fda_03_PMF_full.pdf`
  - `logs/fda_03_de_pca.log`

**Step 4 — `fda_04_de_gsea.py`**
- **Outputs**
  - `figures/fda_04_gsea_targets.pdf`
  - `figures/fda_04_leadingedge_oxdetox_MTX.pdf`
  - `logs/fda_04_de_gsea.log`

---

### PROTACs analysis (steps 1–21)

**Step 1 — `protacs_01_globalvarianceanalysis.py`**
- **Outputs**
  - `figures/protacs_01_global_PCA_scree.pdf`
  - `figures/protacs_01_heatmap.png`
  - `logs/protacs_01_globalvarianceanalysis.log`

**Step 2 — `protacs_02_azmetadata_cluster.py`**
- **Outputs**
  - `figures/protacs_02_class_metaplot.pdf`
  - `figures/protacs_02_degrader_metaplot.pdf`
  - `figures/protacs_02_ligand_metaplot.pdf`
  - `logs/protacs_02_azmetadata_cluster.log`

**Step 3 — `protacs_03_azmetadata_dendrogram.R`**
- **Outputs**
  - `figures/protacs_03_ChemicalSeries_dendrogram_plot.pdf`
  - `logs/protacs_03_azmetadata_dendrogram.log`

**Step 4 — `protacs_04_limma_metadataconstructor.py`**
- **Outputs**
  - `logs/protacs_04_limma_metadataconstructor.log`

**Step 5 — `protacs_05_limma_cluster_0p1.R`**
- **Outputs**
  - `figures/protacs_05_volcanoes_0p1uM.png`
  - `logs/protacs_05_limma_cluster_0p1.log`

**Step 6 — `protacs_06_limma_cluster_1p0.R`**
- **Outputs**
  - `figures/protacs_06_volcanoes_1uM.png`
  - `logs/protacs_06_limma_cluster_1p0.log`

**Step 7 — `protacs_07_limma_cluster_10.R`**
- **Outputs**
  - `figures/protacs_07_volcanoes_10uM.png`
  - `logs/protacs_07_limma_cluster_10.log`

**Step 8 — `protacs_08_limma_drug.R`**
- **Outputs**
  - `logs/protacs_08_limma_drug.log`

**Step 9 — `protacs_09_de_pca.py`**
- **Outputs**
  - `figures/protacs_09_3D_PCA_ChemicalSeries.pdf`
  - `figures/protacs_09_PDF_FDA_HBD_Zcentered.pdf`
  - `logs/protacs_09_de_pca.log`

**Step 10 — `protacs_10_de_gsea.py`**
- **Outputs**
  - `figures/protacs_10_GSEA_ChemicalSeries_1and10uM.pdf`
  - `figures/protacs_10_gocc_mitochondrial_protein_containing_complex_lineplot.pdf`
  - `logs/protacs_10_de_gsea.log`

**Step 11 — `protacs_11_de_stringnetworkenrich.py`**
- **Outputs**
  - `figures/protacs_11_STRINGenrich_PROTACs_OnandOfftarg.pdf`
  - `logs/protacs_11_de_stringnetworkenrich.log`

**Step 12 — `protacs_12_de_regression.py`**
- **Outputs**
  - `figures/protacs_12_regression_models(Inner_mitochondrial_membrane_protein_complex).pdf`
  - `logs/protacs_12_de_regression.log`

**Step 13 — `protacs_13_split.py`**
- **Outputs**
  - `logs/protacs_13_split.log`

**Step 14 — `protacs_14_split_limma.py`**
- **Outputs**
  - `logs/protacs_14_split_limma.log`

**Step 15 — `protacs_15_split_enrich.py`**
- **Outputs**
  - `logs/protacs_15_split_enrich.log`

**Step 16 — `protacs_16_gbdt_train_retrain.py`**
- **Outputs (selection)**
  - `figures/protacs_16_xgb_first-pass_shap_explainer.pdf`
  - `figures/protacs_16_xgb_second-pass_shap_explainer.pdf`
  - `figures/protacs_16_xgb_second-pass_calibration_curve.pdf`
  - `scoring_models/protacs_16_xgb_first-pass-model.json`
  - `scoring_models/protacs_16_xgb_second-pass-model.json`
  - `scoring_models/protacs_16_xgb_second-pass-calibrated-model.pkl`
  - `logs/protacs_16_gbdt_train_retrain.log`


**Step 17 — `protacs_17_gbdt_scores.py`**
- **Outputs**
  - `figures/protacs_17_series15_toxscores.pdf`
  - `logs/protacs_17_gbdt_scores.log`

**Step 18 — `protacs_18_pheatmap_topweights_all.R`**
- **Outputs**
  - `figures/protacs_18_pheatmap_topweights_all.pdf`
  - `logs/protacs_18_pheatmap_topweights_all.log`

**Step 19 — `protacs_19_pheatmap_topweights_analogs.R`**
- **Outputs**
  - `figures/protacs_19_pheatmap_topweights_analogs.pdf`
  - `logs/protacs_19_pheatmap_topweights_analogs.log`

**Step 20 — `protacs_20_scores_v_degradation.py`**
- **Outputs**
  - `figures/protacs_20_bimodal_toxscores.pdf`
  - `figures/protacs_20_barplot_toxscores.pdf`
  - `figures/protacs_20_joint_toxscores.pdf`
  - `figures/protacs_20_scatterplot_toxscores.pdf`
  - `logs/protacs_20_scores_v_degradation.log`

**Step 21 — `protacs_21_stats_wetlab.py`**
- **Outputs**
  - `figures/protacs_21_GI50_dotplot.pdf`
  - `figures/protacs_21_complexII_all.pdf`
  - `figures/protacs_21_compound1_glugal.pdf`
  - `figures/protacs_21_galactose.pdf`
  - `figures/protacs_21_seahorsemax.pdf`
  - `figures/protacs_21_xenograft_ARdeg.pdf`
  - `logs/protacs_21_stats_wetlab.log`
  - `logs/protacs_21_summarystats_wetlab.log`

---

### ML / DeepSearch (step 22)

**Step 22 — `protacs_22_gbdt_deepsearch.py`**
- **Outputs (selection)**
  - `figures/protacs_22_xgb_first-pass_shap_explainer.pdf`
  - `figures/protacs_22_xgb_second-pass_shap_explainer.pdf`
  - `figures/protacs_22_xgb_second-pass_calibration_curve.pdf`
  - `scoring_models/protacs_22_xgb_first-pass-model.json`
  - `scoring_models/protacs_22_xgb_second-pass-model.json`
  - `scoring_models/protacs_22_xgb_second-pass-calibrated-model.pkl`
  - `logs/protacs_22_gbdt_train_retrain.log`