#!/bin/bash

# Advanced script to run biofam. For a basic template see run_basic.sh

###################
## START EDITING ##
###################


# Input files as plain text format
#inFolder="/Users/damienarnol1/Documents/local/pro/PhD/FA/results/results_tmofa/merged/"
## inFiles=( "$inFolder/WT.txt" "$inFolder/KO.txt" )
#inFiles=( "$inFolder/all_data.txt" )
#sampleGroups="$inFolder/z_groups.txt"
inFolder=$1
run_ix=$2
inFiles=( "$1/dfm.txt" )
X_Files=( "$1/positions.txt" )

# Options for the input files
delimiter=" " # delimiter, such as "\t", "" or " "
header_rows=1 # set to 1 if the files contain row names
header_cols=1 # set to 1 if the files contain column names

# Output file path, please use the .hdf5 extension
outFolder=$1
outFile=( "$outFolder/sfa_out_$2.hdf5" )

# Data options
center_features=1   # center the features to zero-mean? (not necessary as long as learnMean=1)
scale_views=1 	    # scale the views to unit variance (not necessary as long as there no massive differences in scale)

# Tell if the multi-view MOFA model is used transposed (1 : Yes, 0 : No)
sample_wise_sparsity=0
sample_wise_noise=0

# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
# likelihoods=( gaussian gaussian )
#likelihoods=( gaussian )
#
## Define view names
## views=( wt ko )
#views=( unique )
likelihoods=( gaussian )

# Define view names
views=( unique )

# Define file with covariates (not implemented yet, please ignore)
# covariatesFile="/tmp/covariates.txt"

# Maximum number of iterations
iter=5 # we recommend to set this to a large enough value (>1000)

# Convergence criterion
# Recommendation: a 'tolerance' of 0.01 is quite strict and can take a bit of time, for initial testing we recommend increasing it to 0.1
tolerance=10.0 # training will stop when the change in the evidence lower bound (deltaELBO) is smaller than 0.01
nostop=0       # if nostop=1 the training will complete all iterations even if the convergence criterion is met

# Define the initial number of factors and how inactive factors are dropped during training.
# The model automatically removes inactive factors during training if they explain a fraction of variance smaller than 'dropR2'
# Recommendation:
# (1) If you remove inactive factors (dropR2>0), then the initial number of factors should be large enough
# (2) If you want to get the most strong drivers of variation then we recommend dropR2 to be at least 0.05 (5%), but if you want to capture more subtle sources of variation you should decrease it to 0.01 (1%) or 0.03 (3%)
factors=10   # initial number of factors
startDrop=1  # initial iteration to start shutting down factors
freqDrop=1 	 # frequency of checking for shutting down factors
dropR2=0.0  # threshold on fraction of variance explained

# Define hyperparameters for the feature-wise spike-and-slab sparsity prior
# learnTheta=( 1 1 1 ) 	# 1 means that sparsity is active whereas 0 means the sparsity is inactivated; each element of the vector corresponds to a view
# initTheta=( 1 1 1 ) 	# initial value of sparsity levels (1 corresponds to a dense model, 0.5 corresponds to factors ); each element of the vector corresponds to a view
startSparsity=1000 		# initial iteration to activate the spike and slab, we recommend this to be significantly larger than 1.

# Learn an intercept term (feature-wise means)?
# Recommendation: always leave it active. If all your views are gaussian you can set this to 0 and center the features, it does not matter.
# But for non-gaussian views we noticed that this is very useful, so set it to 1
learnIntercept=0

# Random seed
seed=0 # if 0, the seed is automatically generated using the current time


####################
## FINISH EDITING ##
####################

# Prepare command
cmd='python ../build_model/entry_point.py
	--delimiter "$delimiter"
	--input-files ${inFiles[@]}
	--output-file $outFile
	--likelihoods ${likelihoods[@]}
	--views ${views[@]}
	--groups ${groups[@]}
	--iter $iter
	--tolerance $tolerance
	--start-sparsity ${startSparsity[@]}
	--factors $factors
	--start-drop $startDrop
	--freq-drop $freqDrop
	--drop-r2 $dropR2
	--seed $seed
	--X_Files ${X_Files[@]}
'
#
#--sampleGroups $sampleGroups
if [[ $header_rows -eq 1 ]]; then cmd="$cmd --header-rows"; fi
if [[ $header_cols -eq 1 ]]; then cmd="$cmd --header-cols"; fi
if [[ $features_in_rows -eq 1 ]]; then cmd="$cmd --features-in-rows"; fi
# if [ -n "$covariatesFile" ]; then cmd="$cmd --covariates-file $covariatesFile"; fi
if [[ $center_features -eq 1 ]]; then cmd="$cmd --center-features"; fi
if [[ $scale_views -eq 1 ]]; then cmd="$cmd --scale-views"; fi
if [[ $nostop -eq 1 ]]; then cmd="$cmd --non-stop"; fi
if [[ $learnIntercept -eq 1 ]]; then cmd="$cmd --learn-intercept"; fi

if [[ $sample_wise_sparsity -eq 1 ]]; then cmd="$cmd --sample-wise-sparsity"; fi
if [[ $sample_wise_noise -eq 1 ]]; then cmd="$cmd --sample-wise-noise"; fi

eval $cmd
