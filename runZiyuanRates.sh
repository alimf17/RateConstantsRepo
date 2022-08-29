#!/bin/bash

name=$1

#for cores in 1 2 3 4 5 
#do

echo $name


day=$(date +"%Y%m%d")


Rscript runInference.R numberOfCores nameOfRun rdsFileArrayOfExperimentalData numberOfStepsForRun OutputFilesLocation OptionalBooleanForOnlyAdjacentTransitions OptionalInitialStateForMCMC

#parallel --bibtex

for letters in A B C D
do 
	sem --id $name$day -j 16 "Rscript runInference.R ${name}epe1$letters TM/TM_epe1trans.rds 2000 > /scratch/sigbio_project_root/sigbio_project7/alimf/bufferRates/${name}epe1$letters 2> /scratch/sigbio_project_root/sigbio_project7/alimf/buffer/${name}epe1${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R ${name}nmt81chp2$letters TM/TM_nmt81chp2trans.rds 2000 > /scratch/sigbio_project_root/sigbio_project7/alimf/bufferRates/${name}nmt81chp2$letters 2> /scratch/sigbio_project_root/sigbio_project7/alimf/buffer/${name}nmt81chp2${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R ${name}nmt81clr3$letters TM/TM_nmt81clr3trans.rds 2000 > /scratch/sigbio_project_root/sigbio_project7/alimf/bufferRates/${name}nmt81clr3$letters 2> /scratch/sigbio_project_root/sigbio_project7/alimf/buffer/${name}nmt81clr3${letters}Err"
	sem --id $name$day -j 16 "Rscript runInference.R ${name}nmt81mit1$letters TM/TM_nmt81mit1trans.rds 2000 > /scratch/sigbio_project_root/sigbio_project7/alimf/bufferRates/${name}nmt81mit1$letters 2> /scratch/sigbio_project_root/sigbio_project7/alimf/buffer/${name}nmt81mit1${letters}Err"

done
 
sem --id $name$day --wait
