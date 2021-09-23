#!/bin/bash

#SBATCH --job-name=lhs   			    # Job name
#SBATCH --mail-type=FAIL           		# Turn on e-mail notification (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=carajb@uw.edu       # To use a different email, update here or input from command line at runtime
#SBATCH --export=all                    # Export environment variables to the batch job session

#SBATCH --account=csde    	 		    # Allocation Definition
#SBATCH --partition=csde
#SBATCH --nodes=1   		 		    # Node resources
#SBATCH --mem=120G   		 		    # Memory per node
#SBATCH --time=00:30:00 	 		    # Wall clock time limit hrs:min:sec

USRNAME=carajb    # SET ME: your username
echo "${USRNAME}"
export USRNAME

DIRPATH=/gscratch/csde/${USRNAME}/HHCoM    # SET ME: path to your HHCoM directory
echo "${DIRPATH}"
export DIRPATH

NSETS=5600    # SET ME: nSets, size of parameter matrix
echo "${NSETS}" 
export NSETS

TCURR=0    # t_curr, always zero for LHS
echo "${TCURR}"
export TCURR

DATE=23Sep21    # SET ME: date identifier of calibration
echo "${DATE}"
export DATE

pwd; hostname; date

echo "Running MATLAB script to get LHS parameter sets."

cd ${DIRPATH}

export TZ="America/Los_Angeles"
module load matlab_2018a
mkdir -p /gscratch/csde/${USRNAME}/$SLURM_JOB_ID
matlab -nodisplay -nosplash -r "calib1_lhs(${NSETS} , ${TCURR} , '${DATE}' , '${USRNAME}')"

rm -r /gscratch/csde/${USRNAME}/$SLURM_JOB_ID 

date
