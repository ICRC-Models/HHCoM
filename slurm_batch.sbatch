#!/bin/bash

#SBATCH --job-name=batch_test_${SETIDX}              # Job name
#SBATCH --workdir=/gscratch/csde/carajb/HHCoM        # Working directory
#SBATCH --mail-type=FAIL           		     # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                    	     # Export environment variables to the batch job session

#SBATCH --account=csde      	                     # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		     # Node resources
#SBATCH --mem=58G   		                     # Memory per node
#SBATCH --time=02:00:00 	                     # Wall clock time limit hrs:min:sec
#SBATCH --qos=MaxJobs10

pwd; hostname; date

echo "Running MATLAB script to estimate negSumLog for subset of parameter sets"
echo "${SETIDX}" 
export SETIDX
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "calib2_sumll4sets(${SETIDX} , ${TCURR} , '${DATE}')"

rm -r /gscratch/csde/carajb/$SLURM_JOB_ID

date
