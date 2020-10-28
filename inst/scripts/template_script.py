
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

# Two formats are allowed for the input data:

# Option 1: a nested list of matrices, where the first index refers to the view and the second index refers to the group.
#           samples are stored in the rows and features are stored in the columns.
# 			Missing values must be filled with NAs, including samples missing an entire view

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

file = "ftp://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz"
data = pd.read_csv(file, sep="\t")

###########################
## Initialise MOFA model ##
###########################


## (1) initialise the entry point ##
ent = entry_point()


## (2) Set data options ##
# - scale_groups: if groups have significantly different ranges, it is good practice to scale each group to unit variance
# - scale_views: if views have significantly different ranges, it is good practice to scale each view to unit variance
ent.set_data_options(
	scale_groups = False, 
	scale_views = False
)


## (3, option 1) Set data using the nested list of matrices format ##
views_names = ["view1","view2"]
groups_names = ["groupA","groupB"]

# samples_names nested list with length NGROUPS. Each entry g is a list with the sample names for the g-th group
# - if not provided, MOFA will fill it with default samples names
# samples_names = (...)

# features_names nested list with length NVIEWS. Each entry m is a list with the features names for the m-th view
# - if not provided, MOFA will fill it with default features names
# features_names = (...)

# ent.set_data_matrix(data, 
# 	views_names = views_names, 
# 	groups_names = groups_names, 
# 	samples_names = samples_names,   
# 	features_names = features_names
# )

# (3, option 2) Set data using a long data frame
ent.set_data_df(data)


## (4) Set model options ##
# - factors: number of factors. Default is K=10
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli"). 
# 		Default is None, and they are infered automatically
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_factors: use ARD prior in the factors? (TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (TRUE if using multiple views)

# Simple (using default values)
ent.set_model_options()

# Advanced (using personalised values)
ent.set_model_options(
	factors = 5, 
	spikeslab_weights = True, 
	ard_factors = True, 
	ard_weights = True
)


## (5) Set training options ##
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow". 
#		For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training.
# 		Default is None, inactive factors are not dropped during training
# - gpu_mode: use GPU mode? this needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed

# Simple (using default values)
ent.set_train_options()

# Advanced (using personalised values)
ent.set_train_options(
	iter = 100, 
	convergence_mode = "fast", 
	startELBO = 1, 
	freqELBO = 1, 
	dropR2 = None, 
	gpu_mode = False, 
	verbose = False, 
	seed = 42
)


## (6, optional) Set stochastic inference options##
# Only recommended with very large sample size (>1e6) and when having access to GPUs
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.25 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)

# Simple (using default values)
# ent.set_stochastic_options()

# Advanced (using personalised values)
# ent.set_stochastic_options(batch_size=0.5, learning_rate=0.75, forgetting_rate=0.5, start_stochastic=10)


####################################
## Build and train the MOFA model ##
####################################

# Build the model 
ent.build()

# Run the model
ent.run()

##################################################################
## (Optional) do dimensionality reduction from the MOFA factors ##
##################################################################

# ent.umap()
# ent.tsne()

####################
## Save the model ##
####################

outfile = "/Users/ricard/data/mofaplus/hdf5/test.hdf5"

# - save_data: logical indicating whether to save the training data in the hdf5 file.
# this is useful for some downstream analysis in R, but it can take a lot of disk space.
ent.save(outfile, save_data=True)

######################################################
## (Optional Extract metrics from the trained model ##
######################################################

# NOTE: downstream analysis is done efficiently with the MOFA2 R package

# Extract factors (per group)
factors = ent.model.nodes["Z"].getExpectation()

# Extract weights (per view)
weights = ent.model.nodes["W"].getExpectation()

# Extract variance explained
r2 = ent.model.calculate_variance_explained()


##################################################
## (Optional Extract metrics from the hdf5 file ##
##################################################

# NOTE: downstream analysis is done efficiently with the MOFA2 R package

import h5py

f = h5py.File(outfile, 'r')

# See the groups that are stored in the hdf5 file
f.keys()

# Extract factors (per group)
f["expectations"]["Z"]["group_0"].value
f["expectations"]["Z"]["group_1"].value

# Extract weights (per view)
f["expectations"]["W"]["view_0"].value
f["expectations"]["W"]["view_1"].value

# Extract variance explained estimates
f["variance_explained"]["r2_per_factor"]
f["variance_explained"]["r2_total"]
