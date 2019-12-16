#!/bin/bash

#SBATCH --job-name=loadUp                               # Job name
#SBATCH --chdir=/gscratch/csde/guiliu/HHCoM        # Working directory
#SBATCH --mail-type=ALL                             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=guiliu@uw.edu
#SBATCH --export=all                                 # Export environment variables to the batch job session

#SBATCH --account=csde                               # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                    # Node resources
#SBATCH --mem=120G                                    # Memory per node
#SBATCH --time=00:10:00                              # Wall clock time limit hrs:min:sec	

pwd; hostname; date

echo "Running MATLAB loadUp script to get calibration targets."
STEPSPERYEAR=6
echo "${STEPSPERYEAR}"
export STEPSPERYEAR

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/guiliu/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "loadUp(${STEPSPERYEAR})"
rm -r /gscratch/csde/guiliu/$SLURM_JOB_ID

date

