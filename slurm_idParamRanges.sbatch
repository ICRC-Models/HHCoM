#!/bin/bash

#BATCH --job-name=idRanges                           # Job name
#SBATCH --mail-type=FAIL                             # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --export=all                                 # Export environment variables to the batch job session

#SBATCH --account=csde                               # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1                                    # Node resources
#SBATCH --mem=120G                                   # Memory per node
#SBATCH --time=00:20:00                              # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB idParamRanges script to get parameter ranges of top fifty best-fitting sets and compare to bounds."
echo "${USRNAME}"
export USRNAME
echo "${DIRPATH}"
git export DIRPATH
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE

cd ${DIRPATH}

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/${USRNAME}/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "idParamRanges(${TCURR} , '${DATE}')"
rm -r /gscratch/csde/${USRNAME}/$SLURM_JOB_ID

date
