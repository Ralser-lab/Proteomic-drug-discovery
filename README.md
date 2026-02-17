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

