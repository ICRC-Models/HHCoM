#!/bin/bash

#SBATCH --job-name=futureSim_mox                    # Job name
#SBATCH --chdir=/gscratch/csde/dpwhite/HHCoM         # Working directory
#SBATCH --mail-type=FAIL                            # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=dpwhite@uw.edu
#SBATCH --export=all                                # Export environment variables to the batch job session

#SBATCH --account=csde                              # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                   # Node resources
#SBATCH --mem=185G                                  # Memory per node
#SBATCH --time=06:00:00                             # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script futureSim for simulation predictions."

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/dpwhite/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "futureSim(0 , [] , [] , [] , [] , 0 , '19May20')"

rm -r /gscratch/csde/dpwhite/$SLURM_JOB_ID

date

