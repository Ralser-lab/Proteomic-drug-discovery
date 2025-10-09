#!/usr/bin/env bash
# Shell script to reproduce all figures, models, and analysis in manuscript
# in chronological order

set -euo pipefail
mkdir -p logs

SCRIPTS=()
STATUS=()
HYPER=${1:-}
HYPER_FILE="configs/${1:-}"

# run(): time stamps and extracts log while running a target script.
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

# Execute scripts sequentially with run().
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

# Print Summary and a produce master log-file. 
echo
echo "===================== SUMMARY ====================="
for i in "${!SCRIPTS[@]}"; do
  printf "%-50s %s\n" "${SCRIPTS[$i]}" "${STATUS[$i]}"
done
echo "==================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] >>> All tasks complete"
touch logs/all_manuscript_steps.done
 
