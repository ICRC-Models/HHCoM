#!/bin/bash

#SBATCH --job-name=calibHistorical                  # Job name
#SBATCH --workdir=/gscratch/csde/carajb/HHCoM       # Working directory
#SBATCH --mail-type=FAIL                            # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                                # Export environment variables to the batch job session

#SBATCH --account=csde                              # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                   # Node resources
#SBATCH --mem=58G                                  # Memory per node
#SBATCH --time=02:00:00                             # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script calibratorAll3 for historical simulation."

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "calibratorAll3"

rm -r /gscratch/csde/carajb/$SLURM_JOB_ID

date

