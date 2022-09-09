#!/bin/bash

name=$1

priorName=$2

#for cores in 1 2 3 4 5 
#do

echo $name


day=$(date +"%Y%m%d")


#Rscript runInference.R numberOfCores nameOfRun rdsFileArrayOfExperimentalData numberOfStepsForRun OutputFilesLocation OptionalBooleanForOnlyAdjacentTransitions OptionalInitialStateForMCMC

#parallel --bibtex

for letters in A B C D
do 
	sem --id $name$day -j 8 "Rscript runInferenceEqMAg20201122.R 4 ${name}nmt81clr3$letters TM/TM_nmt81clr3trans.rds 10000 /anvil/scratch/x-farhat/bufferRates/ F /anvil/scratch/x-farhat/bufferRates/${priorName}nmt81clr3${letters}_savestate.rds > /anvil/scratch/x-farhat/bufferRates/${name}nmt81clr3$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}nmt81clr3${letters}Err"
	sem --id $name$day -j 8 "Rscript runInferenceEqMAg20201122.R 4 ${name}nmt81mit1$letters TM/TM_nmt81mit1trans.rds 10000 /anvil/scratch/x-farhat/bufferRates/  F /anvil/scratch/x-farhat/bufferRates/${priorName}nmt81mit1${letters}_savestate.rds > /anvil/scratch/x-farhat/bufferRates/${name}nmt81mit1$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}nmt81mit1${letters}Err"

done
 
sem --id $name$day --wait
