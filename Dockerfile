# ===============================================================
# Proteomic Drug Discovery Environment
# Base: Bioconductor + R 4.3 + Python 3.11.5 + version-locked libs
# ===============================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL maintainer="Shaon Basu" \
      description="Reproducible R + Python environment for proteomic drug discovery pipeline" \
      version="1.0"

# ------------------------------
# System deps for building Python + SciPy stack + Rust bootstrap
# ------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran \
      libblas-dev liblapack-dev \
      libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev \
      libffi-dev libncursesw5-dev tk-dev libgdbm-dev libnss3-dev liblzma-dev uuid-dev \
      curl ca-certificates wget xz-utils \
  && rm -rf /var/lib/apt/lists/*

# ------------------------------
# Install modern Rust (needed for gseapy build on aarch64)
# ------------------------------
ENV RUSTUP_HOME=/usr/local/rust CARGO_HOME=/usr/local/rust
ENV PATH=/usr/local/rust/bin:$PATH
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable \
 && rustc -V && cargo -V

# ------------------------------
# Build and install Python 3.11.5 from source
# ------------------------------
ENV PY_VER=3.11.5
WORKDIR /tmp
RUN wget -q https://www.python.org/ftp/python/${PY_VER}/Python-${PY_VER}.tgz \
 && tar -xzf Python-${PY_VER}.tgz \
 && cd Python-${PY_VER} \
 && ./configure --enable-optimizations --with-lto --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make altinstall \
 && /usr/local/bin/python3.11 -m ensurepip --upgrade \
 && /usr/local/bin/python3.11 -m pip install --upgrade pip setuptools wheel \
 && ln -sf /usr/local/bin/python3.11 /usr/local/bin/python \
 && ln -sf /usr/local/bin/python3.11 /usr/local/bin/python3 \
 && cd / && rm -rf /tmp/Python-${PY_VER} /tmp/Python-${PY_VER}.tgz

# ------------------------------
# Python packages (version-locked)
# ------------------------------
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
      openpyxl==3.1.2


# ------------------------------
# R helper + CRAN packages (exact versions)
# ------------------------------
RUN R -q -e 'install.packages(c("remotes","BiocManager"), repos="https://cloud.r-project.org", Ncpus=parallel::detectCores())'

RUN R -q -e 'pkgs <- c( \
  "ggplot2@3.5.2","dplyr@1.1.4","tidyr@1.3.1","pheatmap@1.0.13","cowplot@1.2.0", \
  "RColorBrewer@1.1-3","ggnewscale@0.5.2","ape@5.8-1","factoextra@1.0.7","ggfortify@0.4.18"); \
  parse <- function(s) sub("@.*","",s); vers <- function(s) sub(".*@","",s); \
  for (s in pkgs) remotes::install_version(package=parse(s), version=vers(s), upgrade="never", repos="https://cloud.r-project.org")'

# ------------------------------
# Bioconductor (release pinned)
# ------------------------------
RUN R -q -e '\
  BiocManager::install(version = "3.18", ask = FALSE, update = FALSE); \
  BiocManager::install(c( \
    "limma" = "3.58.1", \
    "EnhancedVolcano" = "1.20.0", \
    "ComplexHeatmap" = "2.18.0", \
    "ggtree" = "3.10.1", \
    "ggtreeExtra" = "1.12.0" \
  ), ask = FALSE, update = FALSE)'

# ------------------------------
# Final setup
# ------------------------------
WORKDIR /image
COPY . /image
RUN chmod +x /workflow/CODERUNNER.sh || true

# Helpful: print versions at build time
RUN python --version && pip --version && rustc -V && cargo -V && R --version

CMD ["bash"]
