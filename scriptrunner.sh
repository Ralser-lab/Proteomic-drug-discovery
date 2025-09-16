#!/usr/bin/env bash

DATA_DIR="data"
FIG_DIR="figures"

## TO DO:
# Remove timestamps from all codes 
# Remove directory stuff, as scripts run relative now

# SCRIPT="code/fda_01_globalvarianceanalysis_250304a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/fda_02_limma_drug_250304a.R"
# echo ">>> Running $SCRIPT"
# time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

# SCRIPT="code/fda_03_de_pca_250304a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/fda_04_de_gsea_250304a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/protacs_01_globalvarianceanalysis_250304a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/protacs_02_azmetadata_cluster_250306a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/protacs_03_azmetadata_dendrogram_250306a.R"
# echo ">>> Running $SCRIPT"
# time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

# SCRIPT="code/protacs_04_limma_metadataconstructor_250304a.py"
# echo ">>> Running $SCRIPT"
# time python3 "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR"

# SCRIPT="code/protacs_05_limma_cluster_0p1_250305a.R"
# echo ">>> Running $SCRIPT"
# time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

# SCRIPT="code/protacs_06_limma_cluster_1p0_250305a.R"
# echo ">>> Running $SCRIPT"
# time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

# SCRIPT="code/protacs_07_limma_cluster_10_250305a.R"
# echo ">>> Running $SCRIPT"
# time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

SCRIPT="code/protacs_08_limma_drug_250305a.R"
echo ">>> Running $SCRIPT"
time Rscript "$SCRIPT" --data-dir "$DATA_DIR" --fig-dir "$FIG_DIR" --install-deps false

echo ">>> Done"