.mofapy2_dependencies <- c(
    "h5py==3.1.0",
    "pandas==1.2.1",
	"scikit-learn==0.24.1"
)


#' @importFrom basilisk BasiliskEnvironment
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = "mofapy2==0.6.1")