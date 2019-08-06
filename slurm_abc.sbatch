#!/bin/bash

#SBATCH --job-name=abc                               # Job name
#SBATCH --workdir=/gscratch/csde/carajb/HHCoM        # Working directory
#SBATCH --mail-type=FAIL                             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                                 # Export environment variables to the batch job session

#SBATCH --account=csde                               # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                    # Node resources
#SBATCH --mem=58G                                    # Memory per node
#SBATCH --time=03:00:00                              # Wall clock time limit hrs:min:sec	

pwd; hostname; date

echo "Running MATLAB abc_smc script to get next set of particles."
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE
echo "{NSETS}"
export NSETS

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "abc_smc(${TCURR} , '${DATE}' , ${NSETS})"
rm -r /gscratch/csde/carajb/$SLURM_JOB_ID

date

