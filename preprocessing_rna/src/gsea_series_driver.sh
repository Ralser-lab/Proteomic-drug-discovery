#!/bin/bash -l

set -euo pipefail

cd ..

mkdir -p logs
mkdir -p gsea_output
mkdir -p ranked_lists

cd src/gsea/

echo "Creating ranked lists"

python3 create_ranked_lists.py 

# Loop through all files in ranked_lists
for file in ../../ranked_lists/*.csv; do
    fname=$(basename "$file")

    echo "Running GSEA on: $fname"

    python3 run_gsea.py \
        -c "$fname"

echo "Creating GSEA plots."
python3 gsea_plot.py
    
done