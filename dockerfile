FROM rocker/r-ver:4.1.2

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libfontconfig1-dev libbz2-dev liblzma-dev \
    libharfbuzz-dev libfribidi-dev libcairo2-dev

RUN apt-get install -y llvm-10

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn scanpy

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools', 'tidyverse'))"

# Install spatstat
RUN R --no-echo --no-restore --no-save -e "install.packages(c('spatstat.explore', 'spatstat.geom'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install latest Matrix
RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages(c('remotes', 'devtools'))"

RUN R --no-echo --no-restore --no-save -e "devtools::install_version('Seurat', version = '4.1.0', repos = 'http://cran.us.r-project.org')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

# Install Seurat Data
RUN R --no-echo --no-restore --no-save -e "devtools::install_github('satijalab/seurat-data')"

# Install Azimuth
RUN R --no-echo --no-restore --no-save -e "devtools::install_github('satijalab/azimuth@v0.4.4')"

# Install Harmony
RUN R --no-echo --no-restore --no-save -e "devtools::install_version('harmony', version = '0.1.0', repos = 'http://cran.us.r-project.org')"

# Install Symphony
RUN R --no-echo --no-restore --no-save -e "devtools::install_version('symphony', version = '0.1.0', repos = 'http://cran.us.r-project.org')"

COPY ./ /screfmapping/
COPY ./data/cache_symphony_sct.uwot /home/rstudio/autoimmune_10x/

RUN wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
RUN tar xvf pbmc3k_filtered_gene_bc_matrices.tar.gz

WORKDIR /home/rstudio/autoimmune_10x
# CMD [ "R" ]