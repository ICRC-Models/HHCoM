TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=31Oct19
echo "${DATE}"
export DATE

echo "Running MATLAB script to get matrix size."
sbatch -p csde -A csde slurm_sizeMatrix.sbatch
sleep 300
FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})
echo "${NSETS}" 
export NSETS

echo "Running simulations, first try."
INT=0
for SETIDX in $(seq 1 16 ${NSETS}); do 
export SETIDX
sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs10
INT=$(($INT + 1))
if [ $INT -ge 10 ]; then 
sleep 4800 # pause to be kind to the scheduler
INT=0
fi 
done

echo "Running MATLAB script to identify failed simulations."
sbatch -p csde -A csde slurm_idMissing.sbatch
sleep 300
FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
RERUN=$(<$FILE)
echo "$RERUN"

echo "Re-running failed simulations."
while [ ! -z "$RERUN" ]; do
MISSING=($RERUN)
INT=0
for SETIDX in ${MISSING[@]}; do
export SETIDX
sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs10
INT=$(($INT + 1))
if [ $INT -ge 10 ] || [ $INT -eq ${#MISSING[@]} ]; then 
sleep 4800 # pause to be kind to the scheduler
INT=0
fi 
done

echo "Running MATLAB script to identify failed simulations, again."
sbatch -p csde -A csde slurm_idMissing.sbatch
sleep 300
FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
RERUN=$(<$FILE)
echo "$RERUN"
done

echo "Running MATLAB abc_smc script to get next set of particles."
sbatch -p csde -A csde slurm_abc.sbatch
sleep 21600
 
echo "Running MATLAB idParamRanges script to get ranges of parameters in best-fitting sets."
sbatch -p csde -A csde slurm_idParamRanges.sbatch

