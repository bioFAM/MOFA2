.mofapy2_dependencies <- c(
    "h5py==3.1.0",
    "pandas==1.2.1",
	"scikit-learn==0.24.1",
	"dtw-python==1.1.10"
)

.mofapy2_version <- "0.6.4"

#' @importFrom basilisk BasiliskEnvironment
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = paste0("mofapy2==",.mofapy2_version))