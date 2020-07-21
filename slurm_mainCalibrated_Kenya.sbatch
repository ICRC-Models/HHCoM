#!/bin/bash

#SBATCH --job-name=Kenya_historicalSim_mox                      # Job name
#SBATCH --chdir=/gscratch/csde/guiliu/HPV-HIVacq_model/HHCoM         # Working directory
#SBATCH --mail-type=ALL                                         # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=guiliu@uw.edu
#SBATCH --export=all                                            # Export environment variables to the batch job session

#SBATCH --account=csde                                          # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                               # Node resources
#SBATCH --mem=185G                                              # Memory per node
#SBATCH --time=04:00:00                                         # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script historicalSim for historical simulation."

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/guiliu/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "historicalSim_Kenya(0 , [] , [] , [] , [] , 0 , '10Dec19')"

rm -r /gscratch/csde/guiliu/$SLURM_JOB_ID 

date

