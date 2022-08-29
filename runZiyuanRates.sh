#!/bin/bash

name=$1

#for cores in 1 2 3 4 5 
#do

echo $name


day=$(date +"%Y%m%d")


#Rscript runInference.R numberOfCores nameOfRun rdsFileArrayOfExperimentalData numberOfStepsForRun OutputFilesLocation OptionalBooleanForOnlyAdjacentTransitions OptionalInitialStateForMCMC

#parallel --bibtex

for letters in A B C D
do 
	sem --id $name$day -j 16 "Rscript runInference.R 4 ${name}epe1$letters TM/TM_epe1trans.rds 2000 /anvil/scratch/x-farhat/bufferRates/ > /anvil/scratch/x-farhat/bufferRates/${name}epe1$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}epe1${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R 4 ${name}nmt81chp2$letters TM/TM_nmt81chp2trans.rds 2000 /anvil/scratch/x-farhat/bufferRates/ > /anvil/scratch/x-farhat/bufferRates/${name}nmt81chp2$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}nmt81chp2${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R 4 ${name}nmt81clr3$letters TM/TM_nmt81clr3trans.rds 2000 /anvil/scratch/x-farhat/bufferRates/ > /anvil/scratch/x-farhat/bufferRates/${name}nmt81clr3$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}nmt81clr3${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R 4 ${name}nmt81mit1$letters TM/TM_nmt81mit1trans.rds 2000 /anvil/scratch/x-farhat/bufferRates/ > /anvil/scratch/x-farhat/bufferRates/${name}nmt81mit1$letters 2> /anvil/scratch/x-farhat/bufferRates/${name}nmt81mit1${letters}Err"

done
 
sem --id $name$day --wait
