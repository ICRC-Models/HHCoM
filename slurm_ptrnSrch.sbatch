#!/bin/bash

#SBATCH --job-name=ptrnSrch                          # Job name
#SBATCH --chdir=/gscratch/csde/carajb/HHCoM          # Working directory
#SBATCH --mail-type=FAIL           		     # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                    	     # Export environment variables to the batch job session

#SBATCH --account=csde      	                     # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		     # Node resources
#SBATCH --mem=120G   		                     # Memory per node
#SBATCH --time=90:00:00 	                     # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script to estimate negSumLog for subset of parameter sets"

SETIDX=1
echo "${SETIDX}" 
export SETIDX
TCURR=0
echo "${TCURR}"
export TCURR
DATE=19Nov19
echo "${DATE}"
export DATE

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "lhsPatternSrch_prt2(${SETIDX} , ${TCURR} , '${DATE}')"

rm -r /gscratch/csde/carajb/$SLURM_JOB_ID

date
