#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs

SCRIPTS=()
STATUS=()

# run(): function to execute a script
# Detects the correct interpreter (python3 or Rscript) based on file extension,
# runs the script, logs output (stdout/stderr) with timestamps, and reports status.
run() {
  local runner="$1"
  local script="$2"
  local log="logs/$(basename "${script%.*}").log"
  local start=$(date '+%Y-%m-%d %H:%M:%S')

  echo "[$start] >>> Running $script"
  {
    echo "[$start] >>> START $script"
    time "$runner" "$script"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] >>> FINISH $script"
  } >"$log" 2>&1

  if [[ $? -eq 0 ]]; then
    SCRIPTS+=("$script")
    STATUS+=("OK ($log)")
  else
    SCRIPTS+=("$script")
    STATUS+=("FAIL ($log)")
    exit 1
  fi
}

# Execute scripts with run()
# run python3 code/fda_01_globalvarianceanalysis_250304a.py
# run Rscript code/fda_02_limma_drug_250304a.R
# run python3 code/fda_03_de_pca_250304a.py
# run python3 code/fda_04_de_gsea_250304a.py
# run python3 code/protacs_01_globalvarianceanalysis_250304a.py
# run python3 code/protacs_02_azmetadata_cluster_250306a.py
# run Rscript code/protacs_03_azmetadata_dendrogram_250306a.R
# run python3 code/protacs_04_limma_metadataconstructor_250304a.py
# run Rscript code/protacs_05_limma_cluster_0p1_250305a.R
# run Rscript code/protacs_06_limma_cluster_1p0_250305a.R
# run Rscript code/protacs_07_limma_cluster_10_250305a.R
# run Rscript code/protacs_08_limma_drug_250305a.R
# run python3 code/protacs_09_de_pca_250305a.py
# run python3 code/protacs_10_de_gsea_250305.py
# run python3 code/protacs_11_de_stringnetworkenrich_250305a.py
# run python3 code/protacs_12_de_regression_250305a.py
# run python3 code/protacs_13_split_250311a.py
# run Rscript code/protacs_14_split_limma_250306a.R
# run python3 code/protacs_15_split_enrich_250306a.py
run python3 code/protacs_16_ML_brutalgrid_0p01learn_retrain_250306a.py


# REDUCE heatmap.pdf -> heatmaps.png
# Add grand run time in logfiles 

# Print Summary
echo
echo "===================== SUMMARY ====================="
for i in "${!SCRIPTS[@]}"; do
  printf "%-50s %s\n" "${SCRIPTS[$i]}" "${STATUS[$i]}"
done
echo "==================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] >>> All tasks complete"
