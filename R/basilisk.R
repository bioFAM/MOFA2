library(basilisk)

.mofapy2_dependencies <- c(
 #    "h5py==2.10.0",
 #    "pandas==1.1.2",
 #    "pip==20.2.3",
 #    "requests==2.24.0",
	# "scikit-learn==0.23.2",
 #    "wheel==0.35.1",
 #    "six==1.15.0",
 #    "zipp==3.3.0",
    "scipy==1.5.2",
    "argparse==1.4.0",
	"numpy==1.19.2"
)


#' @importFrom basilisk BasiliskEnvironment
#' 
# setupBasiliskEnv("mofa_env", packages, channels = "conda-forge", pip = NULL)
mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies)

proc <- basiliskStart(mofa_env)
# basiliskRun(proc, function() {})

# mofa_env <- basilisk::BasiliskEnvironment("mofa_env", pkgname="MOFA2",
#                             packages=c(
#                                 "pandas==1.1.2",
#                                 "h5py==2.10.0",
#                                 "scipy==1.5.2",
#                                 "argparse==1.4.0",
#                                 "scikit-learn==0.21.3"),
#                             pip = "mofapy2==0.5.6")
