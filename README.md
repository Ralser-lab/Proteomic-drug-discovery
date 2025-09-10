# HT-MS-DrugRefiner

Drug activity deconvolution using high-throughput screening, ultra-fast proteomics, and machine learning. 

Contains the code to regenerate all figures, training and analytics in 'Proteome guided discovery accurately maps and mitigates toxicity mechanisms of therapeutic androgen receptor degraders'.

## Project Structure 

```
project_root/      
│     
│-- code/                 # main folder with objects, scripts, and notebook
│   ├── device_gradientboostingmachine.py 
│   ├── device_summarystatistics.py    
│   ├── fda_01_globalvarianceanalysis_250324a.py - [33.92s]
│   ├── fda_02_limma_drug_250304a.R - [3.78s]  
│   ├── fda_03_de_pca_250304a.py - [19.35s] 
│   ├── fda_04_de_gsea_250304a.py - [232.38s] 
│   ├── protacs_01_globalvarianceanalysis_250304a.py - [91.12s]
│   ├── protacs_02_azmetadata_cluster_250306a.py - [1.29s]
│   ├── protacs_03_azmetadata_dendrogram_250306a.R - [1.06s]
│   ├── protacs_04_limma_metadataconstructor_250304a.py - [5.94s]
│   ├── protacs_05_limma_cluster_0p1_250305a.R - [16.59s]
│   ├── protacs_06_limma_cluster_1p0_250305a.R - [16.56s]
│   ├── protacs_07_limma_cluster_10_250305a.R - [16.48s]
│   ├── protacs_08_limma_drug_250305a.R - [17.85s]
│   ├── protacs_09_de_pca_250305a.py - [8.87s]
│   ├── protacs_10_de_gsea_250305.py - [509.41s] 
│   ├── protacs_11_de_stringnetworkenrich_250305a.py - [2.71s]
│   ├── protacs_12_de_regression_250305a.py - [5.45s] 
│   ├── protacs_13_split_250311a.py - [2.21s]
│   ├── protacs_14_split_limma_250306a.R - [14.88s] 
│   ├── protacs_15_split_enrich_250306a.py - [1.67s]
│   ├── protacs_16_ML_brutalgrid_0p01learn_retrain_250306a.py - [278.25s]
│   ├── protacs_17_ML_finalmodel_predictions_250307.py - [5.71s]
│   ├── protacs_18_ML_dotplot_250307a.R - [1.38s]
│   ├── protacs_19_ML_dotplot_analog_250307a.R - [1.41s]
│   ├── protacs_20_stats_wetlab_250307a.py - [5.02s]
│   ├── protacs_21_safetyscore_distributions_250522a.ipynb
│ 
│-- data/             
│   ├── README.md         # instructions on sourcedata aquisition
│
│-- figures/              # output folder to regenerate all figs
│   ├── README.md           
│ 
│-- scoring_models/       # output folder for ML models    
│   ├── final_model_250305a.json
│   ├── final_calibrated_model_250305a.pkl
│
│-- requirements.txt      # Python dependencies for virtual / local env
│-- renv.lock             # R environment lockfile
│-- .gitignore  
│-- README.md  

```
## Setup Instructions

### Environments

- Python 3.11.5 (download dependencies in a virtual environment from 'requirements.txt') 

- R 4.3.1 + Bioconductor 3.18 (use renv to restore dependencies from 'renv.lock')

### Installation / Download (~10 seconds)

1. Clone directory in local folder 

   ```sh
   git clone https://github.com/BasuShaon/HT-MS-DrugRefiner.git
   cd HT-MS-DrugRefiner
   pip install -r requirements.txt

2. Download source data from FigShare (URL in data availability statement: 420 MB) 

3. Copy contents into `/data`.

### Running the Code (Total runtime: ~22 minutes)

1. Navigate into `/code`. 

   ```sh
   cd code

2. Execute the scripts in alphanumerical order, starting with:

   ```sh
   Python3 fda_01_globalvarianceanalysis_250304a.py 
   Rscript fda_02_limma_drug_250304a.R
   ...

   *Ignore python modules containing device objects*

3. View scores, intermediate and output files in `/data` & `/figures`
