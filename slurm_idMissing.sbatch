#!/bin/bash

#SBATCH --job-name=idMissing                         # Job name
#SBATCH --chdir=/gscratch/csde/carajb/HHCoM          # Working directory
#SBATCH --mail-type=FAIL                             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                                 # Export environment variables to the batch job session

#SBATCH --account=csde                          # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                    # Node resources
#SBATCH --ntasks-per-node=1                          # CPUs per node
#SBATCH --mem=100G                                   # Memory per node
#SBATCH --time=00:15:00                              # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script to identify failed simulations."
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE
echo "${NSETS}"
export NSETS

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "idMissingSets(${TCURR} , '${DATE}' , ${NSETS})"
rm -r /gscratch/csde/carajb/$SLURM_JOB_ID

date

