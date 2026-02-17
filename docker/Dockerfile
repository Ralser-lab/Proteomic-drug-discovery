# ===============================================================
# Proteomic Drug Discovery Environment
# Base: Bioconductor + R 4.3 + Python 3.11.5 + version-locked libs
# ===============================================================

# 1. Base image (PINNED)
FROM bioconductor/bioconductor_docker:RELEASE_3_18@sha256:8c131494af2c1a4a8ff16656ac486e72fd0d1b135d4c06e54f17bf63cd3b86a3

LABEL maintainer="Shaon Basu" \
      description="Reproducible R + Python environment for proteomic drug discovery pipeline" \
      version="1.1"

# 2. System dependencies for Python + SciPy + Rust bootstrap
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran \
      libblas-dev liblapack-dev \
      libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev \
      libffi-dev libncursesw5-dev tk-dev libgdbm-dev libnss3-dev liblzma-dev uuid-dev \
      curl ca-certificates wget xz-utils \
  && rm -rf /var/lib/apt/lists/*

# 3. Install Rust (pin release version)
ENV RUSTUP_HOME=/usr/local/rust \
    CARGO_HOME=/usr/local/rust \
    PATH=/usr/local/rust/bin:$PATH

RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain 1.80.0 \
 && rustc -V && cargo -V

# 4. Build Python (pin release version)
ENV PY_VER=3.11.5
WORKDIR /tmp

RUN wget -q https://www.python.org/ftp/python/${PY_VER}/Python-${PY_VER}.tgz \
 && tar -xzf Python-${PY_VER}.tgz \
 && cd Python-${PY_VER} \
 && ./configure --enable-optimizations --with-lto --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make altinstall \
 && /usr/local/bin/python3.11 -m ensurepip --upgrade \
 && /usr/local/bin/python3.11 -m pip install --upgrade \
      "pip==24.2" "setuptools==70.0.0" "wheel==0.43.0" \
 && ln -sf /usr/local/bin/python3.11 /usr/local/bin/python \
 && ln -sf /usr/local/bin/python3.11 /usr/local/bin/python3 \
 && cd / && rm -rf /tmp/Python-${PY_VER} /tmp/Python-${PY_VER}.tgz

# 5. Python packages (pin package versions)
RUN python -m pip install --no-cache-dir \
      gseapy==1.0.6 \
      joblib==1.3.2 \
      matplotlib==3.8.1 \
      numpy==1.25.2 \
      pandas==2.1.0 \
      scikit-learn==1.3.0 \
      scipy==1.11.2 \
      seaborn==0.13.2 \
      shap==0.46.0 \
      statsmodels==0.14.0 \
      xgboost==2.0.3 \
      openpyxl==3.1.2 \
      snakemake==9.12.0 \
      optuna==4.5.0

# 6. R helper + CRAN packages (pin package versions)
RUN R -q -e 'install.packages(c("remotes","BiocManager"), repos="https://cloud.r-project.org", Ncpus=parallel::detectCores())'

RUN R -q -e 'pkgs <- c( \
  "ggplot2@3.5.2","dplyr@1.1.4","tidyr@1.3.1","pheatmap@1.0.13","cowplot@1.2.0", \
  "RColorBrewer@1.1-3","ggnewscale@0.5.2","ape@5.8-1","factoextra@1.0.7","ggfortify@0.4.18"); \
  parse <- function(s) sub("@.*","",s); vers <- function(s) sub(".*@","",s); \
  for (s in pkgs) remotes::install_version(package=parse(s), version=vers(s), upgrade="never", repos="https://cloud.r-project.org")'

# 7. Bioconductor (pin release version, install packages)
RUN R -q -e '\
  install.packages(c("BiocManager"), repos="https://cloud.r-project.org"); \
  BiocManager::install(version = "3.18", ask = FALSE, update = FALSE); \
  BiocManager::install(c( \
    "limma", \ 
    "EnhancedVolcano", \
    "ComplexHeatmap", \
    "ggtree", \
    "ggtreeExtra" \
  ), ask = FALSE, update = FALSE)'

# 8. Final setup
WORKDIR /image
COPY . /image
RUN chmod +x /run_all.sh || true

# Print build-time versions for transparency
RUN python --version && pip --version && rustc -V && cargo -V && R --version

CMD ["bash"]
