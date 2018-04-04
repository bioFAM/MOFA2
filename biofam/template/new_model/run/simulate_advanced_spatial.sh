#!/bin/bash

# Advanced script to run biofam. For a basic template see run_basic.sh

###################
## START EDITING ##
###################


factors=5   # initial number of factors
M=1
N=100
D=( 100 )
spatialFact=0.2
noise=1.  # 0.5, 1, 2, 3, 4, 5
sparsity=0.2
outDir='/home/yonatan/PycharmProjects/covar/biofam/biofam/template/new_model/run/spatial'
outFile='/home/yonatan/PycharmProjects/covar/biofam/biofam/template/new_model/run/spatial/simul_spatial.h5'


# Tell if the multi-view MOFA model is used transposed (1 : Yes, 0 : No)
transpose=1

#Choose to sample the positions of the samples to test the covariance prior structure (for any view if transpose = 1)
sample_X=1


# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian )


# Names of views
views=( view_A )



# Random seed
seed=0 # if 0, the seed is automatically generated using the current time

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
	--D ${D[@]}
	--N $N
	--spatialFact $spatialFact
	--noise $noise
	--sparsity $sparsity
	--views $views
'

eval $cmd
