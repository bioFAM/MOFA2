#!/bin/bash

# Advanced script to run biofam. For a basic template see run_basic.sh

###################
## START EDITING ##
###################


factors=10   # initial number of factors
M=2
N=100
D=( 200 200 )
noise=1.  # 0.5, 1, 2, 3, 4, 5
sparsity=0.2
outDir='/home/yonatan/PycharmProjects/biofam2/biofam/biofam/template/new_model/run/test_data/'
outFile='/home/yonatan/PycharmProjects/biofam2/biofam/biofam/template/new_model/run/test_data/simul.h5'


# Tell if the multi-view MOFA model is used transposed (1 : Yes, 0 : No)
transpose_sparsity=1
transpose_noise=1


# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian )


# Names of views
views=( view_A view_B )



# Random seed
seed=0 # if 0, the seed is automatically generated using the current time

####################
## FINISH EDITING ##
####################

# Prepare command
cmd='python3 ../build_model/simulation_entry.py
	--outFile $outFile
	--outDir $outDir
	--likelihoods ${likelihoods[@]}
	--factors $factors
	--M $M
	--D ${D[@]}
	--N $N
	--noise $noise
	--sparsity $sparsity
	--views ${views[@]}
	--seed $seed
'

if [[ $transpose_sparsity -eq 1 ]]; then cmd="$cmd --transpose_sparsity"; fi
if [[ $transpose_noise -eq 1 ]]; then cmd="$cmd --transpose_noise"; fi

eval $cmd
