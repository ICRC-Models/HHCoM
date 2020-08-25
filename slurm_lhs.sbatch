#!/bin/bash

#SBATCH --job-name=lhs   			    # Job name
#SBATCH --chdir=/gscratch/csde/carajb/HHCoM         # Working directory
#SBATCH --mail-type=FAIL           		    # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                                # Export environment variables to the batch job session

#SBATCH --account=csde    	 		    # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		    # Node resources
#SBATCH --mem=120G   		 		    # Memory per node
#SBATCH --time=00:30:00 	 		    # Wall clock time limit hrs:min:sec

NSETS=7644    # nSets
echo "${NSETS}" 
export NSETS

TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=22Apr20Ph2V10
echo "${DATE}"
export DATE

pwd; hostname; date

echo "Running MATLAB script to get LHS parameter sets."

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "calib1_lhs(${NSETS} , ${TCURR} , '${DATE}')"

rm -r /gscratch/csde/carajb/$SLURM_JOB_ID 

date
