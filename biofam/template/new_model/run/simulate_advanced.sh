#!/bin/bash

# Advanced script to run biofam. For a basic template see run_basic.sh

###################
## START EDITING ##
###################

# Output file path, please use the .hdf5 extension
outFile=( "/tmp/test.hdf5" )
outDir='/tmp'

# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian gaussian )

factors=20   # initial number of facotrs
M=3
N=100
D=1000

# Random seed
seed=0 # if 0, the seed is automatically generated using the current time

####################
## FINISH EDITING ##
####################

# Prepare command
cmd='python ../build_model/simulation_entry.py
	--outFile $outFile
	--outDir $outDir
	--likelihoods ${likelihoods[@]}
	--factors $factors
	--seed $seed
	--M $M
	--D $D
	--N $N
'

eval $cmd
