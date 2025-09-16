#!/usr/bin/env bash

DATA_DIR="data"
FIG_DIR="figures"

SCRIPT="code/fda_01_globalvarianceanalysis_250304a.py"
echo ">>> Running $SCRIPT"
time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

SCRIPT="code/fda_02_limma_drug_250304a.R"
echo ">>> Running $SCRIPT"
time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

SCRIPT="code/fda_03_de_pca_250304a.py"
echo ">>> Running $SCRIPT"
time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

SCRIPT="code/fda_04_de_gsea_250304a.py"
echo ">>> Running $SCRIPT"
time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

echo ">>> Done"