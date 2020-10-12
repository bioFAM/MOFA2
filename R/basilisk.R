#' @importFrom basilisk BasiliskEnvironment
mofa_env <- basilisk::BasiliskEnvironment("mofa_env", pkgname="MOFA2",
                            packages=c("numpy==1.19.1",
                                       "pandas==1.1.1",
                                       "h5py==2.10.0",
                                       "scipy==1.5.2",
                                       "argparse==1.4.0",
                                       "scikit-learn==0.23.2",
                                       "mofapy2==0.5.6"
                                       ))