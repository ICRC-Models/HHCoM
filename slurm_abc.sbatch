#!/bin/bash

#SBATCH --job-name=abc                               # Job name
#SBATCH --chdir=/gscratch/csde/willmin/HHCoM          # Working directory
#SBATCH --mail-type=FAIL                             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=willmin@uw.edu
#SBATCH --export=all                                 # Export environment variables to the batch job session

#SBATCH --account=csde                          # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                    # Node resources
#SBATCH --ntasks-per-node=28                         # CPUs per node
#SBATCH --mem=120G                                   # Memory per node
#SBATCH --time=24:00:00                              # Wall clock time limit hrs:min:sec	

pwd; hostname; date

echo "Running MATLAB abc_smc script to get next set of particles."
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE
echo "${NSETS}"
export NSETS

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/willmin/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "abc_smc(${TCURR} , '${DATE}' , ${NSETS})"
rm -r /gscratch/csde/willmin/$SLURM_JOB_ID

date

