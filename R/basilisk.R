# .mofapy2_dependencies <- c(
#     "h5py==3.1.0",
#     "pandas==1.2.1",
#   "scikit-learn==0.24.1",
#   "dtw-python==1.1.10"
# )

.mofapy2_dependencies <- c(
    "python=3.10.5",
    "numpy=1.23.1",
    "scipy=1.8.1",
    "pandas=1.4.3",
    "h5py=3.6.0",
    "scikit-learn=1.1.1",
    "dtw-python=1.2.2"
)

.mofapy2_version <- "0.7.0"

#' @importFrom basilisk BasiliskEnvironment
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = paste0("mofapy2==",.mofapy2_version))