# DIA-MS Preprocessing Pipeline

Data-driven pipeline to build complete proteomes with reduced batch effects from DIA-MS precursor fragment intensities. Output follows a normal probability distribution and feeds directly into the main analysis pipeline.

## Input data

Download input files from Figshare and place them into `./preprocessing_dia/input/`:

**[https://doi.org/10.6084/m9.figshare.30469304](https://doi.org/10.6084/m9.figshare.30469304)**

Two files are required:

| File | Description |
|---|---|
| `SB_PROTAC_prmatrix_240314a.tsv` | Precursor matrix — sample IDs as columns, precursor intensities as rows, indexed by protein and fragment ID |
| `20240314_AF_50-0121_metadata.xlsx` | Metadata — sample IDs as rows, MS batch as a categorical column named `MS.Batch` |

## Execution

Run via Snakemake from the repository root:

```bash
docker run --rm -v "$PWD":/image prot-env snakemake -s snakefiles/preprocessing_dia.smk -j 1
```

Output is written to `./preprocessing_dia/output/` and the processed prmatrix is deposited into `./data/` for use by the main pipeline.

## Pipeline modules

- **QC** — sample and precursor quality filters  
- **Filtration** — missingness thresholds across samples and features  
- **Imputation** — minimum-value imputation on remaining missing entries  
- **Batch correction** — ComBat / limma-based correction on MS batch  
- **maxLFQ** — protein-level intensity summarisation (R wrapper)

## Requirements

- Python 3.11.5  
- R 4.3.1  

All dependencies are bundled in the provided [Dockerfile](../docker/Dockerfile).
