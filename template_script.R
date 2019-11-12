library(MOFA2)

###############
## Load data ##
###############

# the data data can take three formats:
# - List of matrices
# - Long data.frame
# - Seurat object

# data <- (...)

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)


####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)

# model options
# - num_factors: number of factors
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (default TRUE)
# - spikeslab_factors: use spike-slab sparsity prior in the factors? (default FALSE)
# - ard_factors: use ARD prior in the factors? (default TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (always TRUE)
MOFAobject_opts <- get_default_model_options(MOFAobject)

# Training options
# - maxiter: number of iterations
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training
# - gpu_mode: use GPU mode? (needs cupy installed and a functional GPU, see https://cupy.chainer.org/)
# - verbose: verbose mode?
# stochastic: use stochastic inference?
# - seed: random seed
train_opts <- get_default_training_options(MOFAobject)

# Set stochastic inference options
# (Optional, only recommended with very large data sets) 
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.1 to 0.5)

# stochastic_opts <- get_default_stochastic_options(MOFAobject)

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
  # stochastic_options = stochastic_opts
)


#####################
## Train the model ##
#####################

outfile <- "(...)/model.hdf5"
MOFAobject.trained <- run_mofa(MOFAobject, outfile)
