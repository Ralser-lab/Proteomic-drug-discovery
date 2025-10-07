# Proteomic-drug-discovery

Repository containing scripts to regenerate all figures, training and analytics in 
**'Proteome-guided discovery accurately maps and mitigates toxicity mechanisms of therapeutic androgen receptor degraders'.**

## Installation

1. Clone this repository and navigate to root directory:
```bash
git clone https://github.com/BasuShaon/Proteomic-drug-discovery.git
cd /Proteomic-drug-discovery
```
2. Download data [here](https://figshare.com/s/6d164fd50adfdb9a68d7) and copy-paste it into `/data`
   
## Reproduce Environment 

1. Install [Docker](https://www.docker.com/get-started)
2. Build environment image with docker (copy-paste in root directory):
```bash
   docker build -t proteomediscovery-env .
```
### Details:

**Python 3.11.5** (gseapy==1.0.6, joblib==1.3.2, matplotlib==3.8.1, numpy==1.25.2, openpyxl==3.1.2, pandas==2.1.0, scikit-learn==1.3.0, scipy==1.11.2, seaborn==0.13.2, shap==0.46.0, statsmodels==0.14.0, xgboost==2.0.3) 

**R 4.3.1** (ggplot2==3.5.2, dplyr==1.1.4, tidyr==1.3.1, pheatmap==1.0.13, cowplot==1.2.0, RColorBrewer==1.1-3, ggnewscale==0.5.2, ape==5.8-1, factoextra==1.0.7, ggfortify==0.4.18) 

**Bioconductor 3.18** (limma==3.58.1, EnhancedVolcano=1.20.0, ComplexHeatmap==2.18.0, ggtree==3.10.1, ggtreeExtra=1.12.0) 

**Tested Hardware** MacBook Pro (M2 MAX CPU, 32 GB RAM, macOS Ventura 13.3) 

## Execution (~20 minute runtime)

Run & load `CODERUNNER.sh` in a docker container within the root directory (copy-paste this):
```bash
   docker run -it --rm -v "$PWD":/image proteomediscovery-env bash CODERUNNER.sh HYPER.json
```

To adjust hyperparameter grid for the xgboost toxicity scoring workflow, edit `HYPER.json` with desired search space, save, and run CODERUNNER.sh with docker in command-line as described above. 

Generated outputs (models, figures, logfiles) save into `/scoring_models`, `/figures`, and `/logs` respectively, and they map as follows:

### Mapping index

```
project_root/      
│     
│-- code/            
│   ├── device_gradientboostingmachine.py 
│   ├── device_summarystatistics.py  
│
│   ├── fda_01_globalvarianceanalysis
│   │   ├── figures/fda_01_dispersion_kde.pdf
│   │   ├── figures/fda_01_global_PCA.pdf
│   │   ├── figures/fda_01_global_PCA_scree.pdf
│   │   ├── figures/fda_01_heatmap.pdf
│   │   ├── figures/fda_01_heatmap.png
│   │   ├── logs/fda_01_globalvarianceanalysis.log
│
│   ├── fda_02_limma_drug
│   │   ├── figures/fda_02_volcanoes_top3.png
│   │   ├── logs/fda_02_limma_drug.log
│
│   ├── fda_03_de_pca
│   │   ├── figures/fda_03_2D_PCA_tmatrix.pdf
│   │   ├── figures/fda_03_3D_PCA_tmatrix.pdf
│   │   ├── figures/fda_03_PDF_all.pdf
│   │   ├── figures/fda_03_PDF_stats.pdf
│   │   ├── figures/fda_03_PMF_conditioned.pdf
│   │   ├── figures/fda_03_PMF_full.pdf
│   │   ├── logs/fda_03_de_pca.log
│
│   ├── fda_04_de_gsea
│   │   ├── figures/fda_04_gsea_targets.pdf
│   │   ├── figures/fda_04_leadingedge_oxdetox_MTX.pdf
│   │   ├── logs/fda_04_de_gsea.log
│
│   ├── protacs_01_globalvarianceanalysis
│   │   ├── figures/protacs_01_dispersion_kde.pdf
│   │   ├── figures/protacs_01_global_PC1_GSEA.pdf
│   │   ├── figures/protacs_01_global_PC2_GSEA.pdf
│   │   ├── figures/protacs_01_global_PC3_GSEA.pdf
│   │   ├── figures/protacs_01_global_PCA_GSEA.pdf
│   │   ├── figures/protacs_01_global_PCA_MSBatch.pdf
│   │   ├── figures/protacs_01_global_PCA_concentration.pdf
│   │   ├── figures/protacs_01_global_PCA_scree.pdf
│   │   ├── figures/protacs_01_heatmap.png
│   │   ├── logs/protacs_01_globalvarianceanalysis.log
│
│   ├── protacs_02_azmetadata_cluster
│   │   ├── figures/protacs_02_class_metaplot.pdf
│   │   ├── figures/protacs_02_degrader_metaplot.pdf
│   │   ├── figures/protacs_02_ligand_metaplot.pdf
│   │   ├── logs/protacs_02_azmetadata_cluster.log
│
│   ├── protacs_03_azmetadata_dendrogram
│   │   ├── figures/protacs_03_ChemicalSeries_dendrogram_plot.pdf
│   │   ├── logs/protacs_03_azmetadata_dendrogram.log
│
│   ├── protacs_04_limma_metadataconstructor
│   │   ├── logs/protacs_04_limma_metadataconstructor.log
│
│   ├── protacs_05_limma_cluster_0p1
│   │   ├── figures/protacs_05_volcanoes_0p1uM.png
│   │   ├── logs/protacs_05_limma_cluster_0p1.log
│
│   ├── protacs_06_limma_cluster_1p0
│   │   ├── figures/protacs_06_volcanoes_1uM.png
│   │   ├── logs/protacs_06_limma_cluster_1p0.log
│
│   ├── protacs_07_limma_cluster_10
│   │   ├── figures/protacs_07_volcanoes_10uM.png
│   │   ├── logs/protacs_07_limma_cluster_10.log
│
│   ├── protacs_08_limma_drug
│   │   ├── logs/protacs_08_limma_drug.log
│
│   ├── protacs_09_de_pca
│   │   ├── figures/protacs_09_3D_PCA_ChemicalSeries.pdf
│   │   ├── figures/protacs_09_PDF_FDA_HBD_Zcentered.pdf
│   │   ├── figures/protacs_09_PDF_violinbox_FDA_HBD.pdf
│   │   ├── figures/protacs_09_PMF_ARHBD.pdf
│   │   ├── figures/protacs_09_PMF_ontargetHBD.pdf
│   │   ├── logs/protacs_09_de_pca.log
│
│   ├── protacs_10_de_gsea
│   │   ├── figures/protacs_10_GSEA_ChemicalSeries_1and10uM.pdf
│   │   ├── figures/protacs_10_PMF_GSEA_PROTACs_all.pdf
│   │   ├── figures/protacs_10_gocc_mitochondrial_protein_containing_complex_lineplot.pdf
│   │   ├── logs/protacs_10_de_gsea.log
│
│   ├── protacs_11_de_stringnetworkenrich
│   │   ├── figures/protacs_11_STRINGenrich_PROTACs_OnandOfftarg.pdf
│   │   ├── logs/protacs_11_de_stringnetworkenrich.log
│
│   ├── protacs_12_de_regression
│   │   ├── figures/protacs_12_barplot_IC50.pdf
│   │   ├── figures/protacs_12_regression_NDUFA5.pdf
│   │   ├── figures/protacs_12_regression_NDUFA5_pred.pdf
│   │   ├── figures/protacs_12_regression_models(Inner_mitochondrial_membrane_protein_complex).pdf
│   │   ├── figures/protacs_12_regression_p_values.pdf
│   │   ├── logs/protacs_12_de_regression.log
│
│   ├── protacs_13_split
│   │   ├── logs/protacs_13_split.log
│
│   ├── protacs_14_split_limma
│   │   ├── logs/protacs_14_split_limma.log
│
│   ├── protacs_15_split_enrich
│   │   ├── logs/protacs_15_split_enrich.log
│
│   ├── protacs_16_gbdt_train_retrain
│   │   ├── figures/protacs_16_xgb_first-pass_classif_report.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_confusion_all_0.51.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_confusion_test_0.51.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_confusion_train_0.51.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_f1thresh.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_pr.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_roc.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_shap_explainer.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_shap_interactions.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_shap_interactions_heatmap.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_top10_features_by_gain.pdf
│   │   ├── figures/protacs_16_xgb_first-pass_top10_features_by_weight.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_classif_report.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_confusion_all_0.67.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_confusion_test_0.67.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_confusion_train_0.67.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_f1thresh.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_pr.pdf
│   │   ├── figures/protacs_16_xgb_second-pass-calibrated_roc.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_calibration_curve.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_classif_report.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_confusion_all_0.55.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_confusion_test_0.55.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_confusion_train_0.55.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_f1thresh.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_pr.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_roc.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_shap_explainer.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_shap_interactions.pdf
│   │   ├── figures/protacs_16_xgb_second-pass_shap_interactions_heatmap.pdf
│   │   ├── scoring_models/protacs_16_xgb_first-pass-final_metrics.csv
│   │   ├── scoring_models/protacs_16_xgb_first-pass-model.json
│   │   ├── scoring_models/protacs_16_xgb_first-pass-search_space.csv
│   │   ├── scoring_models/protacs_16_xgb_second-pass-calibrated-model.pkl
│   │   ├── scoring_models/protacs_16_xgb_second-pass-final_metrics.csv
│   │   ├── scoring_models/protacs_16_xgb_second-pass-model.json
│   │   ├── scoring_models/protacs_16_xgb_second-pass-search_space.csv
│   │   ├── logs/protacs_16_gbdt_train_retrain.log
│
│   ├── protacs_17_gbdt_scores
│   │   ├── figures/protacs_17_analogues_signature_PCA.pdf
│   │   ├── figures/protacs_17_analogues_signature_barplot.pdf
│   │   ├── figures/protacs_17_series15_toxscores.pdf
│   │   ├── logs/protacs_17_gbdt_scores.log
│
│   ├── protacs_18_pheatmap_topweights_all
│   │   ├── figures/protacs_18_pheatmap_topweights_all.pdf
│   │   ├── logs/protacs_18_pheatmap_topweights_all.log
│
│   ├── protacs_19_pheatmap_topweights_analogs
│   │   ├── figures/protacs_19_pheatmap_topweights_analogs.pdf
│   │   ├── logs/protacs_19_pheatmap_topweights_analogs.log
│
│   ├── protacs_20_scores_v_degradation
│   │   ├── figures/protacs_20_barplot_toxscores.pdf
│   │   ├── figures/protacs_20_bimodal_toxscores.pdf
│   │   ├── figures/protacs_20_joint_toxscores.pdf
│   │   ├── figures/protacs_20_scatterplot_toxscores.pdf
│   │   ├── logs/protacs_20_scores_v_degradation.log
│
│   ├── protacs_21_stats_wetlab
│   │   ├── figures/protacs_21_GI50_dotplot.pdf
│   │   ├── figures/protacs_21_complexII_all.pdf
│   │   ├── figures/protacs_21_compound1_glugal.pdf
│   │   ├── figures/protacs_21_galactose.pdf
│   │   ├── figures/protacs_21_seahorsemax.pdf
│   │   ├── figures/protacs_21_xenograft_ARdeg.pdf
│   │   ├── logs/protacs_21_stats_wetlab.log
│   │   ├── logs/protacs_21_summarystats_wetlab.log
│
│   ├── protacs_22_gbdt_deepsearch
│   │   ├── figures/protacs_22_xgb_first-pass_classif_report.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_confusion_all_0.51.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_confusion_test_0.51.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_confusion_train_0.51.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_f1thresh.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_pr.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_roc.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_shap_explainer.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_shap_interactions.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_shap_interactions_heatmap.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_top10_features_by_gain.pdf
│   │   ├── figures/protacs_22_xgb_first-pass_top10_features_by_weight.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_classif_report.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_confusion_all_0.67.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_confusion_test_0.67.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_confusion_train_0.67.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_f1thresh.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_pr.pdf
│   │   ├── figures/protacs_22_xgb_second-pass-calibrated_roc.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_calibration_curve.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_classif_report.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_confusion_all_0.55.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_confusion_test_0.55.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_confusion_train_0.55.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_f1thresh.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_pr.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_roc.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_shap_explainer.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_shap_interactions.pdf
│   │   ├── figures/protacs_22_xgb_second-pass_shap_interactions_heatmap.pdf
│   │   ├── scoring_models/protacs_22_xgb_first-pass-final_metrics.csv
│   │   ├── scoring_models/protacs_22_xgb_first-pass-model.json
│   │   ├── scoring_models/protacs_22_xgb_first-pass-search_space.csv
│   │   ├── scoring_models/protacs_22_xgb_second-pass-calibrated-model.pkl
│   │   ├── scoring_models/protacs_22_xgb_second-pass-final_metrics.csv
│   │   ├── scoring_models/protacs_22_xgb_second-pass-model.json
│   │   ├── scoring_models/protacs_22_xgb_second-pass-search_space.csv
│   │   ├── logs/protacs_22_gbdt_train_retrain.log

