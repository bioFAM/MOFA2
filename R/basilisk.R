.mofapy2_dependencies <- c(
    # "python==3.9",
    "h5py==2.10.0",
    "pandas==1.1.3",
    "pip==20.2.3",
    "requests==2.23.0",
    "wheel==0.34.2",
    "six==1.14.0",
    "zipp==0.3.3",
	"scikit-learn==0.23.2",
    "scipy==1.5.2",
    # "argparse==1.4.0", (in Python by default)
	"numpy==1.19.2",
    if (basilisk.utils::isMacOSX()) "nomkl==3.0"
)


#' @importFrom basilisk BasiliskEnvironment
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = "mofapy2==0.5.6")