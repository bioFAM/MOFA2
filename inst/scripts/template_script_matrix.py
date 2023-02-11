
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

# The data format is a nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be explicitly filled using NAs, including samples missing an entire view

datadir = "/Users/ricard/data/mofaplus/test"
views = ["0","1"]
groups = ["0","1"]
data = [None]*len(views)
for m in range(len(views)):
    data[m] = [None]*len(groups)
    for g in range(len(groups)):
        datafile = "%s/%s_%s.txt.gz" % (datadir, views[m], groups[g])
        data[m][g] = pd.read_csv(datafile, header=None, sep=' ')

###########################
## Initialise MOFA model ##
###########################

## (1) initialise the entry point
ent = entry_point()


## (2) Set data options
# - scale_views: if views have very different ranges, one can to scale each view to unit variance
ent.set_data_options(
	scale_views = False
)


## (3) Define names
views_names = ["view1","view2"]
# groups_names = ["groupA","groupB"]

# samples_names nested list with length n_groups. Each entry g is a list with the sample names for the g-th group
# - if not provided, MOFA will fill it with default samples names
samples_names = (...)

# features_names nested list with length NVIEWS. Each entry m is a list with the features names for the m-th view
# - if not provided, MOFA will fill it with default features names
features_names = (...)


## (4) Set data matrix
ent.set_data_matrix(data, 
	views_names = views_names, 
	groups_names = groups_names, 
	samples_names = samples_names,   
	features_names = features_names
)


## (5) Set model options
# - factors: number of factors. Default is 15
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli"). Default and recommended is "gaussian"
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_weights: use automatic relevance determination prior in the weights? (TRUE if using multiple views)

# using default values
ent.set_model_options()

# using personalised values
ent.set_model_options(
	factors = 5, 
	spikeslab_weights = True, 
	ard_weights = True
)

## (5) Set training options ##
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow". Fast mode is usually good enough.
# - dropR2: minimum variance explained criteria to drop factors while training. Default is None, inactive factors are not dropped during training
# - gpu_mode: use GPU mode? this functionality needs cupy installed and a functional GPU, see https://biofam.github.io/MOFA2/gpu_training.html
# - seed: random seed

# using default values
ent.set_train_options()

# using personalised values
ent.set_train_options(
	iter = 100, 
	convergence_mode = "fast", 
	dropR2 = None, 
	gpu_mode = False, 
	seed = 42
)

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

outfile = "/Users/ricard/data/mofaplus/hdf5/test.hdf5"

# - save_data: logical indicating whether to save the training data in the hdf5 file.
# this is useful for some downstream analysis in R, but it can take a lot of disk space.
ent.save(outfile, save_data=True)

#########################
## Downstream analysis ##
#########################

# Check the mofax package for the downstream analysis in Python: https://github.com/bioFAM/mofax
# Check the MOFA2 R package for the downstream analysis in R: https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html
# All tutorials: https://biofam.github.io/MOFA2/tutorials.html

# Extract factor values (a list with one matrix per sample group)
factors = ent.model.nodes["Z"].getExpectation()

# Extract weights  (a list with one matrix per view)
weights = ent.model.nodes["W"].getExpectation()

# Extract variance explained values
r2 = ent.model.calculate_variance_explained()

# Interact directly with the hdf5 file
import h5py
f = h5py.File(outfile, 'r')
f.keys()

# Extract factors
f["expectations"]["Z"]["group_0"].value
f["expectations"]["Z"]["group_1"].value

# Extract weights
f["expectations"]["W"]["view_0"].value
f["expectations"]["W"]["view_1"].value

# Extract variance explained estimates
f["variance_explained"]["r2_per_factor"]
f["variance_explained"]["r2_total"]
