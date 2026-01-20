#!/bin/bash -l

#SBATCH -J limma
#SBATCH -o "/home/kgkd182/projects/SafetyOmicsDili/preprocessing/logs/limma_%j.out"
#SBATCH -e "/home/kgkd182/projects/SafetyOmicsDili/preprocessing/logs/limma_%j.err"
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mail-user=kgkd182
#SBATCH --mail-type=END,FAIL

set -e

# Load R
module load R/4.3.2-foss-2021a   

# Go to where the script lives 
cd /home/kgkd182/projects/SafetyOmicsDili/preprocessing

# Set data directory
DATA_DIR="data_v2"
mkdir -p ../${DATA_DIR}

echo 'combat...'
Rscript 01_combat.R --data-dir ${DATA_DIR}

echo 'filter...'
Rscript 02_filter_transform.R --data-dir ${DATA_DIR}

echo 'limma...'
Rscript 03_limma_all.R --norm none --data-dir ${DATA_DIR}
Rscript 03_limma_all.R --norm tmm --data-dir ${DATA_DIR}
Rscript 03_limma_all.R --norm voom --data-dir ${DATA_DIR}

echo 'extract...'
Rscript 04_matrix_all.R --norm none --data-dir ${DATA_DIR}
Rscript 04_matrix_all.R --norm tmm --data-dir ${DATA_DIR}
Rscript 04_matrix_all.R --norm voom --data-dir ${DATA_DIR}

echo 'gene2symbol...'
Rscript 05_gene2symbol.R --norm none --data-dir ${DATA_DIR}
Rscript 05_gene2symbol.R --norm tmm --data-dir ${DATA_DIR}
Rscript 05_gene2symbol.R --norm voom --data-dir ${DATA_DIR}

echo "everything saved in ../${DATA_DIR}"
