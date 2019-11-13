
######################################################
## Template script to train a MOFA+ model in Python ##
######################################################

from mofapy2.run.entry_point import entry_point
import pandas as pd
import io
import requests # to download the online data

###############
## Load data ##
###############

# Multiple formats are allowed for the input data:

# Option 1: a nested list of matrices, where the first index refers to the view and the second index refers to the group.
#           samples are stored in the rows and features are stored in the columns.
#           Importantly, all views for a given group G must have the same samples in the rows.
#           If there is any sample that is missing a particular view, the column needs to be filled with NAs

# datadir = "/Users/ricard/data/mofaplus/test"
# views = ["0","1"]
# groups = ["0","1"]
# data = [None]*len(views)
# for m in range(len(views)):
#     data[m] = [None]*len(groups)
#     for g in range(len(groups)):
#         datafile = "%s/%s_%s.txt.gz" % (datadir, views[m], groups[g])
#         data[m][g] = pd.read_csv(datafile, header=None, sep=' ')

# Option 2: a data.frame with columns ["sample","feature","view","group","value"]
#           In this case there is no need to have missing values in the data.frame,
#           they will be automatically filled in when creating the corresponding matrices

file = "ftp://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/mofa2/getting_started/data.txt.gz" # Simulated data
data = pd.read_csv(file, sep="\t")

###########################
## Initialise MOFA model ##
###########################

# initialise the entry point
ent = entry_point()

# Set data options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli")
# - scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance
# - scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance
ent.set_data_options(likelihoods=["gaussian","gaussian"], scale_groups=False, scale_views=False)

# Set data (option 1, nested list of matrices) 
# ent.set_data_matrix(data)

# Set data (option 2, data.frame)
ent.set_data_df(data)

# Set model options
# - factors: number of factors
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli")
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_factors: use ARD prior in the factors? (set to TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (always TRUE)
ent.set_model_options(factors=5, likelihoods=["gaussian","gaussian"], spikeslab_weights=True, ard_factors=True, ard_weights=True)

# Set training options
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - elbofreq: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training
# - gpu_mode: use GPU mode? (needs cupy installed and a functional GPU, see https://cupy.chainer.org/)
# - verbose: verbose mode?
# - seed: random seed
ent.set_train_options(iter=1000, convergence_mode="medium", startELBO=1, elbofreq=1, dropR2=None, gpu_mode=False, verbose=False, seed=42)

# (Optional) Set stochastic inference options
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.1 to 0.5)

####################################
## Build and train the MOFA model ##
####################################

# Build the model 
ent.build()

# Run the model
ent.run()

####################
## Save the model ##
####################

outfile="/Users/ricard/data/mofaplus/hdf5/test.hdf5"
ent.save(outfile)
