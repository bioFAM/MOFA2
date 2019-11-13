library(MOFA2)
library(data.table)

###############
## Load data ##
###############

# Multiple formats are allowed for the input data:

# -- Option 1 --
# Option 1: 
# a nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Importantly, all views for a given group G must have the same samples in the rows.
# If there is any sample that is missing a particular view, the column needs to be filled with NAs

# datadir = "/Users/ricard/data/mofaplus/test"
# views = ["0","1"]
# groups = ["0","1"]
# data = [None]*len(views)
# for m in range(len(views)):
#     data[m] = [None]*len(groups)
#     for g in range(len(groups)):
#         datafile = "%s/%s_%s.txt.gz" % (datadir, views[m], groups[g])
#         data[m][g] = pd.read_csv(datafile, header=None, sep=' ')

# -- Option 2 --
# data.frame with columns ["sample","feature","view","group","value"]
# In this case there is no need to have missing values in the data.frame,
# they will be automatically filled in when creating the corresponding matrices

file = "ftp://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/mofa2/getting_started/data.txt.gz" # Simulated data
data = fread(file)

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
# - scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance (default FALSE)
# - scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance
# - views: view names
# - groups: group names
data_opts <- get_default_data_options(MOFAobject)

# model options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli")
# - num_factors: number of factors
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (default TRUE)
# - spikeslab_factors: use spike-slab sparsity prior in the factors? (default FALSE)
# - ard_factors: use ARD prior in the factors? (default TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (always TRUE)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

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
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42

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

outfile <- paste0(getwd(),"/model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)
