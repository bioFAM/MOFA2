FROM r-base:4.0.2

WORKDIR /mofa2
ADD . /mofa2

RUN apt-get update && apt-get install -f && apt-get install -y python3 python3-setuptools python3-dev
RUN apt-get install -y libcurl4-openssl-dev 
RUN apt-get install -y libcairo2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libxt-dev libharfbuzz-dev libfribidi-dev
RUN python3 setup.py install

# Install bioconductor dependencies
RUN R --vanilla -e "\
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cran.r-project.org'); \
  sapply(c('rhdf5', 'dplyr', 'tidyr', 'reshape2', 'pheatmap', 'corrplot', \
           'ggplot2', 'ggbeeswarm', 'scales', 'GGally', 'doParallel', 'RColorBrewer', \
           'cowplot', 'ggrepel', 'foreach', 'reticulate', 'HDF5Array', 'DelayedArray', \
           'ggpubr', 'forcats', 'Rtsne', 'uwot', \
           'systemfonts', 'ragg', 'Cairo', 'ggrastr'), \ 
         BiocManager::install)"
RUN R CMD INSTALL --build .

CMD []
