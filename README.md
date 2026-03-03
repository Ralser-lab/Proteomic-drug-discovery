# Proteomic-drug-discovery

Repository containing scripts to regenerate all figures, training and analytics in 
**'Proteome-guided discovery accurately maps and mitigates toxicity mechanisms of therapeutic androgen receptor degraders'.**

## Reproduce manuscript findings

1. Clone this repository and navigate to root directory.
2. Download [data](https://doi.org/10.6084/m9.figshare.28578113) and place it into ./data.
3. Install [docker](https://www.docker.com/get-started) and add to command line interface (CLI) $PATH.
4. Build a image of the computing environment using the provided [dockerfile](./docker/Dockerfile).
5. Reproduce analysis by executing `run_all.sh` in a docker container (copy-paste in CLI):
```bash
   git clone https://github.com/BasuShaon/Proteomic-drug-discovery
   cd Proteomic-drug-discovery
   docker build -t prot-env -f docker/Dockerfile . 
   docker run --rm -v "$PWD":/image prot-env bash run_all.sh
```

## Machine learning workflow: Hyperparameter search 

To replicate or adjust the hyperparameters used for ML training: 
1. Edit [./configs/HYPER.json](./configs/HYPER.json).
2. Run `deepsearch.smk` in a docker container (copy-paste in CLI):
```bash
   docker run --rm -v "$PWD":/image prot-env snakemake -s snakefiles/deepsearch.smk -j 12
```

## Pre-processing workflow: Data-independent-aquisition mass spectrometry

To reproduce pre-processing pipeline on DIA-NN prmatrix:
1. Download [input files](https://doi.org/10.6084/m9.figshare.30469304.v1) and place it into ./preprocessing-dia/input.
2. Run `preprocessing-dia.smk` in a docker container (copy-paste in CLI):
```bash
   docker run --rm -v "$PWD":/image prot-env snakemake -s snakefiles/preprocessing-dia.smk -j 12
```

## Pre-processing workflow: Next generation sequencing

To reproduce pre-processing pipeline on Illumina sequencing outputs:
1. Download [data]() and place it into ./preprocessing-ngs/data`.
2. Run `preprocessing-ngs.smk` in a docker container (copy-paste in CLI):
```bash
   docker run --rm -v "$PWD":/image prot-env snakemake -s snakefiles/preprocessing-ngs.smk -j 12
```

## Environment

**Hardware** (MacBook Pro, M2 MAX 12-core CPU, 32 GB RAM, macOS Ventura 13.3) 

The following [software environment](./docker/Dockerfile) was used for development:

**Python 3.11.5** (gseapy==1.0.6, joblib==1.3.2, matplotlib==3.8.1, numpy==1.25.2, openpyxl==3.1.2, optuna==4.5.0, pandas==2.1.0, scikit-learn==1.3.0, scipy==1.11.2, seaborn==0.13.2, shap==0.46.0, snakemake==9.12.0, statsmodels==0.14.0, xgboost==2.0.3) 

**R 4.3.1** (ggplot2==3.5.2, dplyr==1.1.4, tidyr==1.3.1, pheatmap==1.0.13, cowplot==1.2.0, RColorBrewer==1.1-3, ggnewscale==0.5.2, ape==5.8-1, factoextra==1.0.7, ggfortify==0.4.18) 

**Bioconductor 3.18** (limma==3.58.1, EnhancedVolcano==1.20.0, ComplexHeatmap==2.18.0, ggtree==3.10.1, ggtreeExtra==1.12.0) 

## Repository Structure
```bash
├── LICENSE
├── README.md
├── code
│   ├── device_gradientboostingmachine.py
│   ├── device_summarystatistics.py
│   ├── device_supportfunctions.py
│   ├── fda_01_globalvarianceanalysis.py
│   ├── fda_02_limma_drug.R
│   ├── fda_03_de_pca.py
│   ├── fda_04_de_gsea.py
│   ├── protacs_01_globalvarianceanalysis.py
│   ├── protacs_02_azmetadata_cluster.py
│   ├── protacs_03_azmetadata_dendrogram.R
│   ├── protacs_04_limma_metadataconstructor.py
│   ├── protacs_05_limma_cluster_0p1.R
│   ├── protacs_06_limma_cluster_1p0.R
│   ├── protacs_07_limma_cluster_10.R
│   ├── protacs_08_limma_drug.R
│   ├── protacs_09_de_pca.py
│   ├── protacs_10_de_gsea.py
│   ├── protacs_11_de_stringnetworkenrich.py
│   ├── protacs_12_de_regression.py
│   ├── protacs_13_split.py
│   ├── protacs_14_split_limma.R
│   ├── protacs_15_split_enrich.py
│   ├── protacs_16_gbdt_train_retrain.py
│   ├── protacs_17_gbdt_train_wide.py
│   ├── protacs_18_gbdt_scores.py
│   ├── protacs_19_pheatmap_topweights_all.R
│   ├── protacs_20_pheatmap_topweights_analogs.R
│   ├── protacs_21_scores_v_degradation.py
│   ├── protacs_22_stats_wetlab.py
│   ├── protacs_23_rna_gsea.py
│   └── protacs_24_rna_gex.py
├── configs
│   ├── cfg_gridcv.json
│   └── cfg_optuna_wide.json
├── docker
│   └── Dockerfile
├── figures
│   ├── fda_01_dispersion_kde.pdf
│   ├── fda_01_raw_heatmap.png
│   ├── fda_02_volcanoes.png
│   ├── fda_03_de_kde.pdf
│   ├── fda_03_de_pca.pdf
│   ├── fda_04_gsea_targets.pdf
│   ├── fda_04_leadingedge_oxdetox_MTX.pdf
│   ├── protacs_01_dispersion_kde.pdf
│   ├── protacs_01_raw_heatmap.png
│   ├── protacs_01_raw_pc1_gsea.pdf
│   ├── protacs_01_raw_pc2_gsea.pdf
│   ├── protacs_01_raw_pc3_gsea.pdf
│   ├── protacs_01_raw_pca_MSBatch.pdf
│   ├── protacs_01_raw_pca_concentration.pdf
│   ├── protacs_01_raw_pca_gsea.pdf
│   ├── protacs_01_raw_pca_scree.pdf
│   ├── protacs_02_class_metaplot.pdf
│   ├── protacs_02_degrader_metaplot.pdf
│   ├── protacs_02_ligand_metaplot.pdf
│   ├── protacs_03_chemicalseries_dendrogram_plot.pdf
│   ├── protacs_05_volcanoes_0p1uM.png
│   ├── protacs_06_volcanoes_1uM.png
│   ├── protacs_07_volcanoes_10uM.png
│   ├── protacs_09_de_kde_fda_v_hbd_boxplot.png
│   ├── protacs_09_de_kde_fda_v_hbd_zcentered.png
│   ├── protacs_09_de_pca_chemicalseries.pdf
│   ├── protacs_10_de_gsea_chemicalseries_1and10uM.pdf
│   ├── protacs_10_de_gsea_pmf.pdf
│   ├── protacs_10_gocc_mitochondrial_protein_containing_complex_lineplot.pdf
│   ├── protacs_11_de_stringdb_onoff_targ.pdf
│   ├── protacs_12_barplot_IC50.pdf
│   ├── protacs_12_barplot_IC50_mean.pdf
│   ├── protacs_12_barplot_logIC50_mean.pdf
│   ├── protacs_12_regression_NDUFA5.pdf
│   ├── protacs_12_regression_models.pdf
│   ├── protacs_16_xgb_first-pass_pr.pdf
│   ├── protacs_16_xgb_first-pass_shap_explainer.pdf
│   ├── protacs_16_xgb_first-pass_shap_interactions.pdf
│   ├── protacs_16_xgb_first-pass_shap_interactions_heatmap.pdf
│   ├── protacs_16_xgb_first-pass_top10_features_by_gain.pdf
│   ├── protacs_16_xgb_first-pass_top10_features_by_weight.pdf
│   ├── protacs_16_xgb_second-pass-calibrated_pr.pdf
│   ├── protacs_16_xgb_second-pass_calibration_curve.pdf
│   ├── protacs_16_xgb_second-pass_pr.pdf
│   ├── protacs_16_xgb_second-pass_shap_explainer.pdf
│   ├── protacs_16_xgb_second-pass_shap_interactions.pdf
│   ├── protacs_16_xgb_second-pass_shap_interactions_heatmap.pdf
│   ├── protacs_17_xgb_first-pass_pr.pdf
│   ├── protacs_17_xgb_first-pass_shap_explainer.pdf
│   ├── protacs_17_xgb_first-pass_shap_interactions.pdf
│   ├── protacs_17_xgb_first-pass_shap_interactions_heatmap.pdf
│   ├── protacs_17_xgb_first-pass_top10_features_by_gain.pdf
│   ├── protacs_17_xgb_first-pass_top10_features_by_weight.pdf
│   ├── protacs_17_xgb_second-pass-calibrated_pr.pdf
│   ├── protacs_17_xgb_second-pass_calibration_curve.pdf
│   ├── protacs_17_xgb_second-pass_pr.pdf
│   ├── protacs_17_xgb_second-pass_shap_explainer.pdf
│   ├── protacs_17_xgb_second-pass_shap_interactions.pdf
│   ├── protacs_17_xgb_second-pass_shap_interactions_heatmap.pdf
│   ├── protacs_18_analogues_signature_PCA.pdf
│   ├── protacs_18_analogues_signature_barplot.pdf
│   ├── protacs_18_series15_toxscores.pdf
│   ├── protacs_19_pheatmap_topweights_all.pdf
│   ├── protacs_20_pheatmap_topweights_analogs.pdf
│   ├── protacs_21_barplot_toxscores.pdf
│   ├── protacs_21_joint_KDE_toxscores.pdf
│   ├── protacs_21_scatterplot_toxscores_v_IC50.pdf
│   ├── protacs_21_violin_toxscores.pdf
│   ├── protacs_22_GI50_dotplot.pdf
│   ├── protacs_22_cmpd1_glugal.pdf
│   ├── protacs_22_cmpd2_glugal.pdf
│   ├── protacs_22_cmpd3_glugal.pdf
│   ├── protacs_22_complexII_all.pdf
│   ├── protacs_22_galactose.pdf
│   ├── protacs_22_seahorsemax.pdf
│   ├── protacs_22_xenograft_ARdeg.pdf
│   ├── protacs_23_rna_C4_2_gsea_dot_plot.pdf
│   ├── protacs_23_rna_LNCaP_gsea_dot_plot.pdf
│   ├── protacs_24_rna_C4_2_24h_HALLMARK_ANDROGEN_RESPONSE.pdf
│   └── protacs_24_rna_LNCaP_24h_HALLMARK_ANDROGEN_RESPONSE.pdf
├── pyproject.toml
├── run_all.sh
├── scoring_models
│   ├── protacs_16_xgb_first-pass-final_metrics.csv
│   ├── protacs_16_xgb_first-pass-model.json
│   ├── protacs_16_xgb_first-pass-search_space.csv
│   ├── protacs_16_xgb_second-pass-calibrated-model.pkl
│   ├── protacs_16_xgb_second-pass-final_metrics.csv
│   ├── protacs_16_xgb_second-pass-model.json
│   ├── protacs_16_xgb_second-pass-search_space.csv
│   ├── protacs_17_xgb_first-pass-final_metrics.csv
│   ├── protacs_17_xgb_first-pass-model.json
│   ├── protacs_17_xgb_first-pass-search_space.csv
│   ├── protacs_17_xgb_second-pass-calibrated-model.pkl
│   ├── protacs_17_xgb_second-pass-final_metrics.csv
│   ├── protacs_17_xgb_second-pass-model.json
│   └── protacs_17_xgb_second-pass-search_space.csv
└── snakefiles
    ├── ml_train_wide.smk
    ├── preprocessing_dia.smk
    └── preprocessing_rna.smk
```
