#!/bin/bash

#SBATCH --job-name=ptrnSrch_${SETIDX}              # Job name
#SBATCH --chdir=/gscratch/csde/guiliu/HHCoM        # Working directory
#SBATCH --mail-type=ALL           		     # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=guiliu@uw.edu
#SBATCH --export=all                    	     # Export environment variables to the batch job session

#SBATCH --account=csde      	                     # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		     # Node resources
#SBATCH --mem=120G   		                     # Memory per node
#SBATCH --time=6:00:00 	                     # Wall clock time limit hrs:min:sec
#SBATCH --qos=MaxJobs4

pwd; hostname; date

echo "Running MATLAB script to estimate negSumLog for subset of parameter sets"
SETIDX=1
echo "${SETIDX}" 
export SETIDX
TCURR=0
echo "${TCURR}"
export TCURR
DATE=19NOV19
echo "${DATE}"
export DATE

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/guiliu/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "lhsPatternSrch_prt2(${SETIDX} , ${TCURR} , '${DATE}')"

rm -r /gscratch/csde/guiliu/$SLURM_JOB_ID

date
