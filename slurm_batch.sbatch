#!/bin/bash

#SBATCH --job-name=batch                             # Job name
#SBATCH --chdir=/gscratch/csde/guiliu/HPV-HIVacq_model/HHCoM          # Working directory
#SBATCH --mail-type=ALL           		             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=guiliu@uw.edu
#SBATCH --export=all                    	         # Export environment variables to the batch job session

#SBATCH --account=csde    	                         # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		                 # Node resources
#SBATCH --mem=185G   		                         # Memory per node
#SBATCH --time=12:00:00 	                         # Wall clock time limit hrs:min:sec

pwd; hostname; date  #--qos=MaxJobs4

echo "Running MATLAB script to estimate negSumLog for subset of parameter sets"
TCURR=0
DATE=14DEC20
SETIDX=1
echo "${SETIDX}" 
export SETIDX
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/guiliu/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "calib2_sumll4sets(${SETIDX} , ${TCURR} , '${DATE}')"

rm -r /gscratch/csde/guiliu/$SLURM_JOB_ID

date
