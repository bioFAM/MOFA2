FROM r-base:3.5.2

WORKDIR /biofam
ADD . /biofam

RUN apt-get update && apt-get install -y python3 python3-setuptools
RUN python3 setup.py install

# Install bioconductor dependencies
RUN R --vanilla -e "\
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
  sapply(c('rhdf5', 'dplyr', 'tidyr', 'reshape2', 'pheatmap', 'corrplot', \
           'ggplot2', 'ggbeeswarm', 'scales', 'GGally', 'doParallel', 'RColorBrewer', \
           'cowplot', 'ggrepel', 'foreach', 'reticulate', 'HDF5Array', 'DelayedArray'), \ 
         BiocManager::install)"
RUN R CMD INSTALL --build BioFAMtools

CMD []
