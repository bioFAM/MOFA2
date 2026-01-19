# .mofapy2_dependencies <- c(
#     "h5py==3.1.0",
#     "pandas==1.2.1",
#     "scikit-learn==0.24.1",
#     "dtw-python==1.1.10"
# )

.mofapy2_dependencies <- c(
    "python=3.12.12",
    "numpy=1.26.4",
    "scipy=1.12.0",
    "pandas=2.2.1",
    "h5py=3.10.0",
    "scikit-learn=1.4.0",
    "dtw-python=1.3.1"
)

.mofapy2_version <- "0.7.3"

#' @importFrom basilisk BasiliskEnvironment
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = paste0("mofapy2==",.mofapy2_version))