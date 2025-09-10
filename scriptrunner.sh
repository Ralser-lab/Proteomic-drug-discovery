#!/usr/bin/env bash 

echo ">>> Running fda_01_globalvarianceanalysis_250304a.py"
time python3 code/fda_01_globalvarianceanalysis_250304a.py

echo ">>> Running fda_02_limma_drug_250304a.R"
time Rscript code/fda_02_limma_drug_250304a.R

echo ">>> Running fda_03_de_pca_250304a.py"
time python3 code/fda_03_de_pca_250304a.py

echo ">>> Running fda_04_de_gsea_250304a.py"
time python3 code/fda_04_de_gsea_250304a.py

echo ">>> Done"