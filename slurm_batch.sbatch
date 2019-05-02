#!/bin/bash

#SBATCH --workdir=/gscratch/csde/carajb/myjobdir            # Working directory
#SBATCH --mail-type=ALL           							# Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                    					# Export environment variables to the batch job session

#SBATCH --account=csde      	 # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 # Node resources
#SBATCH --mem=58G   		 # Memory per node
#SBATCH --time=01:30:00 	 # Wall clock time limit hrs:min:sec
pwd; hostname; date

echo “Running MATLAB script”

Run script to get sets of sampled parameters
Loop through sets
-feed 16 indices at a time into script that creates a job and calls calibratorAll3 and runs those 16 in parallel, saves results
-note: set # in matlab to 16
export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/xyz/abc/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r mainCalibrated (PARAM1)
date