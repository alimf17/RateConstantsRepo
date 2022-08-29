#!/bin/bash


#SBATCH --account=sigbio_project7
#SBATCH --partition=sigbio
#SBATCH --job-name=ZiyuanRun1  # job name
#SBATCH --output=%x_%j.log             # stdout fname (will be <job-name>_<job_id??>.log) Not sure about job_id here to be honest. It's numbers... just numbers.
#SBATCH --error=%x_%j.err              # stderr fname
#SBATCH --mail-type=BEGIN,END,FAIL     # will get email at job start, finish, and if it fails
#SBATCH --mail-user=alimf@umich.edu # email (please don't use mine...)
#SBATCH --nodes=1                      # Total number of nodes
#SBATCH --ntasks-per-node=64           # Total # of mpi tasks
#SBATCH --time=72:00:00                 # Run time (hh:mm:ss)

# The application(s) to execute along with its input arguments and options:
module load R
module load parallel


bash runZiyuanRates.sh $SLURM_JOB_NAME

