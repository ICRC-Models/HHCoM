#!/bin/bash

#SBATCH --job-name=ptrnSrch                          # Job name
#SBATCH --mail-type=FAIL           		     # Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --export=all                    	     # Export environment variables to the batch job session

#SBATCH --account=csde      	                     # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		     # Node resources
#SBATCH --mem=120G   		                     # Memory per node
#SBATCH --time=90:00:00 	                     # Wall clock time limit hrs:min:sec

pwd; hostname; date

echo "Running MATLAB script to estimate negSumLog for subset of parameter sets"

echo "${USRNAME}"
export USRNAME
echo "${DIRPATH}"
export DIRPATH
echo "${SETIDX}" 
export SETIDX
echo "${TCURR}"
export TCURR
echo "${DATE}"
export DATE

cd ${DIRPATH}

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/${USRNAME}/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "lhsPatternSrch_prt2(${SETIDX} , ${TCURR} , '${DATE}' , '${USRNAME}')"
rm -r /gscratch/csde/${USRNAME}/$SLURM_JOB_ID

date
