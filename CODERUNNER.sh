
#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs

SCRIPTS=()
STATUS=()
HYPER_FILE="${1:-}"

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

# run_hyper(): for scripts that need the HYPER.json arg
run_hyper() {
  local runner="$1"
  local script="$2"
  local log="logs/$(basename "${script%.*}").log"
  local start=$(date '+%Y-%m-%d %H:%M:%S')

  echo "[$start] >>> Running $script with $HYPER_FILE"
  {
    echo "[$start] >>> START $script with $HYPER_FILE"
    time "$runner" "$script" "$HYPER_FILE"
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
run python3 code/fda_01_globalvarianceanalysis.py
run Rscript code/fda_02_limma_drug.R
run python3 code/fda_03_de_pca.py
run python3 code/fda_04_de_gsea.py
run python3 code/protacs_01_globalvarianceanalysis.py
run python3 code/protacs_02_azmetadata_cluster.py
run Rscript code/protacs_03_azmetadata_dendrogram.R
run python3 code/protacs_04_limma_metadataconstructor.py
run Rscript code/protacs_05_limma_cluster_0p1.R
run Rscript code/protacs_06_limma_cluster_1p0.R
run Rscript code/protacs_07_limma_cluster_10.R
run Rscript code/protacs_08_limma_drug.R
run python3 code/protacs_09_de_pca.py
run python3 code/protacs_10_de_gsea.py
run python3 code/protacs_11_de_stringnetworkenrich.py
run python3 code/protacs_12_de_regression.py
run python3 code/protacs_13_split.py
run Rscript code/protacs_14_split_limma.R
run python3 code/protacs_15_split_enrich.py
run python3 code/protacs_16_gbdt_train_retrain.py
run python3 code/protacs_17_gbdt_scores.py
run Rscript code/protacs_18_pheatmap_topweights_all.R
run Rscript code/protacs_19_pheatmap_topweights_analogs.R
run python3 code/protacs_20_scores_v_degradation.py
run python3 code/protacs_21_stats_wetlab.py
run_hyper python3 code/protacs_22_gbdt_deepsearch.py

# Print Summary
echo
echo "===================== SUMMARY ====================="
for i in "${!SCRIPTS[@]}"; do
  printf "%-50s %s\n" "${SCRIPTS[$i]}" "${STATUS[$i]}"
done
echo "==================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] >>> All tasks complete"
 
