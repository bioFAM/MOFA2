#!/bin/bash

# Advanced script to run biofam. For a basic template see run_basic.sh

###################
## START EDITING ##
###################


factors=10   # initial number of factors
M=2
N=100
D=200
spatialFact=0.
noise=1.  # 0.5, 1, 2, 3, 4, 5
sparsity=0.2
outDir='/home/yonatan/PycharmProjects/biofam/biofam/biofam/template/new_model/run/train_non_spatial'
outFile='/home/yonatan/PycharmProjects/biofam/biofam/biofam/template/new_model/run/train_non_spatial/simul_spatial.h5'


# Tell if the multi-view MOFA model is used transposed (1 : Yes, 0 : No)
transpose=0

#Choose to sample the positions of the samples to test the covariance prior structure (for any view if transpose = 1)
sample_X=0


# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian )


# Names of views
views=( view_A view_B )



# Random seed
# seed=2018 # if 0, the seed is automatically generated using the current time

####################
## FINISH EDITING ##
####################

# Prepare command
cmd='python ../build_model/simulation_entry.py
    --transpose $transpose
    --sample_X $sample_X
	--outFile $outFile
	--outDir $outDir
	--likelihoods ${likelihoods[@]}
	--factors $factors
	--M $M
	--D $D
	--N $N
	--spatialFact $spatialFact
	--noise $noise
	--sparsity $sparsity
	--views ${views[@]}
'

eval $cmd
