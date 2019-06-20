#!/bin/bash

#SBATCH --job-name=lhs_test   			    # Job name
#SBATCH --workdir=/gscratch/csde/carajb/HHCoM       # Working directory
#SBATCH --mail-type=FAIL           		    # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu
#SBATCH --export=all                                # Export environment variables to the batch job session

#SBATCH --account=csde      	 		    # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		    # Node resources
#SBATCH --mem=58G   		 		    # Memory per node
#SBATCH --time=00:30:00 	 		    # Wall clock time limit hrs:min:sec

#SBATCH –-output=./SlurmOut/%j.out      # Standard output
#SBATCH –-error=./SlurmOut/%j.err       # Standard error

NSETS=10000    # nSets
echo "${NSETS}" 
export NSETS


pwd; hostname; date

echo "Running MATLAB script to get LHS parameter sets."

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/carajb/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "sensitivityAnalysis3(${NSETS})"

rm -r /gscratch/csde/carajb/$SLURM_JOB_ID 

date
