#!/bin/bash

# Basic script to run MOFA. For a more advanced template see run_advanced.sh

###################
## START EDITING ##
###################

# Input files as plain text format
inFolder="/Users/ricard/mofa/mofa/run/test_data"
inFiles=( "$inFolder/500_0.txt" "$inFolder/500_1.txt" "$inFolder/500_2.txt")

# Options for the input files
delimiter=" " # Delimiter, such as "\t", "" or " "
header_rows=0 # Set to 1 if the files contain row names
header_cols=0 # Set to 1 if the files contain column names

# Output file path, please use the .hdf5 extension
outFile=( "/Users/ricard/MOFA/MOFA/test/test.hdf5" )

# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian gaussian )

# Define view names
views=( A B C )

# Maximum number of iterations
# Recommendation: set this to large enough value (>1000)
iterations=5000

# Define the initial number of factors and how inactive factors are dropped during training.
# The model automatically removes inactive factors during training if they explain a fraction of variance smaller than 'dropR2'
# Recommendation: 
# (1) If you remove inactive factors (dropR2>0), then the initial number of factors should be large enough
# (2) If you want to get the most strong drivers of variation then we recommend dropR2 to be at least 0.05 (5%), but if you want to capture more subtle sources of variation you should decrease it to 0.01 (1%) or 0.03 (3%)
factors=25
dropR2=0.05

####################
## FINISH EDITING ##
####################

# Prepare command
cmd='mofa
	--delimiter "$delimiter"
	--inFiles ${inFiles[@]}
	--outFile $outFile
	--likelihoods ${likelihoods[@]}
	--views ${views[@]}
	--iter $iterations
	--factors $factors
	--dropR2 $dropR2
	--learnIntercept
'
if [[ $header_rows -eq 1 ]]; then cmd="$cmd --header_rows"; fi
if [[ $header_cols -eq 1 ]]; then cmd="$cmd --header_cols"; fi

# Run!
eval $cmd

