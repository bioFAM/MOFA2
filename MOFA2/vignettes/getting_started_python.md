---
title: "MOFA+: training a model in Python"
author: "Ricard Argelaguet"
date: "`r Sys.Date()`"
---

This notebook contains a detailed tutorial on how to train MOFA using Python. If you arrived here, that means that you followed the Python path. Good job, but make sure to respect those that prefer R. They are human beings.  

A concise template script can be found [here](https://github.com/bioFAM/MOFA2/blob/master/template_script.py)

## 1) Load libraries

```{python }
from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

# initialise the entry point
ent = entry_point()
```

## 2) Load data

To create a MOFA+ object you need to specify four dimensions: samples (cells), features, view(s) and group(s). MOFA objects can be created from a wide range of input formats:

### 2.1) pandas data.frame format
A pandas data.frame with columns ["sample","group","feature","view","value"].  
This is personally my favourite format, as it summarises all omics/groups in a single data structure. Also, there is no need to add rows that correspond to missing data.

For example:
```
sample  group   feature value   view
sample1  groupA    gene1   2.8044  RNA
sample1  groupA    gene3    2.2069  RNA
sample2  groupB    gene2    0.1454  RNA
sample2  groupB    gene3     2.7021  RNA
sample2  groupB    promoter1    3.8618  Methylation
sample3  groupB    promoter2    3.2545  Methylation
sample3  groupB    promoter3   1.5014  Methylation
```

Here we load a simulated data set with the following dimensions:
```python
D = [1000,1000] # Number of features per view
M = len(D)      # Number of views
K = 5           # Number of factors
N = [100,100]   # Number of samples per group
G = len(N)      # Number of groups

data_dt = pd.read_csv("/Users/ricard/data/mofaplus/test/data.txt.gz", sep="\t")
data_dt.head()
```

### 2.2) List of matrices
A nested list of numpy arrays, where the first index refers to the view and the second index refers to the group. Samples are stored in the rows and features are stored in the columns. All views for a given group G must have the same samples in the rows. If there is any sample that is missing a 
particular view, the column needs to be filled with NAs.

Loading the same data above in matrix format:
```{python }
data_prefix = "/Users/ricard/data/mofaplus/test"
data_mat = [[None for g in range(G)] for m in range(M)]
for m in range(M):
    for g in range(G):
        data_mat[m][g] = np.loadtxt("%s/%d_%d.txt.gz" % (data_prefix,m,g), delimiter="\t")
```


### 2.3 Define data options

- **likelihoods**: likelihood per view (options are "gaussian", "poisson", "bernoulli")
- **scale_groups**: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is False
- **scale_views**: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is False
```{python }
lik = ["gaussian", "gaussian"]

ent.set_data_options(
    likelihoods = lik, 
    scale_groups = False, 
    scale_views = False
)
```

### 2.4) Add the data to the model

This has to be run after defining the data options
```{python }
# option 1: data.frame format (slower)
ent.set_data_df(data_dt)

# option 2: nested matrix format (faster)
ent.set_data_matrix(data_mat)

# AnnData format
# (...)
```

## 3) Set model options

- **factors**: number of factors
- **likelihods**: same as in data_opts
- **spikeslab_weights**: use spike-slab sparsity prior in the weights? default is TRUE
- **ard_factors**: use ARD prior in the factors? Default is TRUE if using multiple groups.
- **ard_weights**: use ARD prior in the weights? Default is TRUE. 

Only change the default model options if you are familiar with the underlying mathematical model!

```{python }
ent.set_model_options(
    factors = 10, 
    likelihoods = lik, 
    spikeslab_weights = True, 
    ard_factors = True,
    ard_weights = True
)
```

## 4) Set training options

- **iter**: number of iterations. Default is 1000.
- **convergence_mode**: "fast", "medium", "slow". For exploration, the fast mode is good enough.
- **startELBO**: initial iteration to compute the ELBO (the objective function used to assess convergence)
- **elbofreq**: frequency of computations of the ELBO (the objective function used to assess convergence)
- **dropR2**: minimum variance explained criteria to drop factors while training
gpu_mode: use GPU mode? (needs cupy installed and a functional GPU, see https://cupy.chainer.org/)
- **verbose**: verbose mode?
- **seed**: random seed

```{python }
ent.set_train_options(
    iter = 1000, 
    convergence_mode = "fast", 
    startELBO = 1, 
    elbofreq = 1, 
    dropR2 = 0.001, 
    gpu_mode = True, 
    verbose = False, 
    seed = 1
)
```

# 5) (Optional)  stochastic inference options

If the number of samples is very large (at the order of >1e4), you may want to try the stochastic inference scheme. If combined with GPUs, it makes inference significantly faster. However, it requires some additional hyperparameters that in some data sets may need to be optimised (vignette in preparation).

- **batch_size**: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
- **learning_rate**: learning rate (we recommend values from 0.25 to 0.75)
- **forgetting_rate**: forgetting rate (we recommend values from 0.1 to 0.5)

```{python }
# We do not want to use stochastic inference for this data set
<!-- ent.set_stochastic_options(
  batch_size = 0.5,
  learning_rate = 0.5, 
  forgetting_rate = 0.25
) -->
```

## 6) Build and train the MOFA object 

After training, the model will be saved as an hdf5 file
```{python }
ent.build()

ent.run()

# Save the output
outfile = "/Users/ricard/data/mofaplus/hdf5/test.hdf5"
ent.save(outfile)
```

If everything is successful, you should observe an output analogous to the following:
```

######################################
## Training the model with seed 1 ##
######################################

Iteration 1: time=0.03, ELBO=-52650.68, deltaELBO=837116.802 (94.082647669%), Factors=10

(...)

Iteration 9: time=0.04, ELBO=-50114.43, deltaELBO=23.907 (0.002686924%), Factors=10

#######################
## Training finished ##
#######################

Saving model in /Users/ricard/data/mofa2/hdf5/model.hdf5...
```


## Downstream analysis

This finishes the tutorial on how to train a MOFA object from python. To continue with the downstream analysis you have to switch to R. Please, follow [this tutorial](XXX)
