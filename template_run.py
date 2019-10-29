from mofa2py.run.entry_point import entry_point
import pandas as pd

###############
## Load data ##
###############

# Multiple formats are allowed for the input data:

# Option 1: a nested list of matrices, where the first index refers to the view and the second index refers to the group.
#           samples are stored in the rows and features are stored in the columns.
#           Importantly, all views for a given group G must the same number of samples.
#           If there is any sample that is missing a particular view, the column needs to be filled with NAs

# Option 2: a data.frame with columns ["sample","feature","view","group","value"]
#           In this case there is no need to have missing values in the data.frame,
#           they will be automatically filled in when creating the corresponding matrices

datadir = "/Users/ricard/data/mofaplus/test"

views = ["0","1"]
sample_groups = ["0","1"]
data = [None]*len(views)
for m in range(len(views)):
    data[m] = [None]*len(sample_groups)
    for g in range(len(sample_groups)):
        datafile = "%s/test_%s_%s.txt" % (datadir, views[m], sample_groups[g])
        data[m][g] = pd.read_csv(datafile, header=None, sep=' ')

# Dewfine likelihoods
lik = ["gaussian"]*len(data)

###########################
## Initialise MOFA model ##
###########################

# initialise the entry point
ent = entry_point()

# Set data options
# - likelihoods:
# - scale_groups:
# - scale_views
ent.set_data_options(likelihoods=lik, scale_groups=False, scale_views=False)

# Set data (option 1, nested list of matrices)
ent.set_data_matrix(data)

# Set data (option 2, data.frame)
# ent.set_data_df(data)

# Set model options
# - factors
# - likelihods
# - spikeslab_factors
# - spikeslab_weights
# - ard_factors
# - ard_weights
ent.set_model_options(factors=5, likelihoods=lik, spikeslab_factors=False, spikeslab_weights=True, ard_factors=True, ard_weights=True)

# Set training options
# - iter
# - convergence_mode
# - startELBO
# - elbofreq
# - dropR2
# - gpu_mode
# - verbose
# - seed
ent.set_train_options(iter=10, convergence_mode="fast", startELBO=1, elbofreq=1, dropR2=None, gpu_mode=True, verbose=False, seed=1)

# (Optional) Set stochastic vb options
# - batch_size
# - forgetting_rate
# - learning_rate
# ent.set_stochastic_options(batch_size=.5, learning_rate=0.75, forgetting_rate=0.5)

####################################
## Build and train the MOFA model ##
####################################

ent.build()
ent.run()
ent.impute()

####################
## Save the model ##
####################

outfile="/Users/ricard/data/mofaplus/hdf5/test.hdf5"
ent.save(outfile)
