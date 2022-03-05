USRNAME=carajb    # SET ME: your username
echo "${USRNAME}"
export USRNAME

DIRPATH=/gscratch/csde/${USRNAME}/HHCoM    # SET ME: path to your HHCoM directory
echo "${DIRPATH}"
export DIRPATH

TCURR=6    # t_curr, last iteration of calibration
echo "${TCURR}"
export TCURR

DATE=22Apr20Ph2V11    # date identifier of calibration
echo "${DATE}"
export DATE

FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})    # size of parameter matrix
echo "${NSETS}" 
export NSETS

echo "Running historical simulations."
SETIDX=1
export SETIDX
sbatch -p csde -A csde --mail-user=${USRNAME}@uw.edu slurm_runMultHistSims.sbatch    # note: values passed from command line have precedence over values defined in job script

