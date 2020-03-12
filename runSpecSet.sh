TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=24Feb20
echo "${DATE}"
export DATE

echo "Get matrix size."
FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})
echo "${NSETS}" 
export NSETS

echo "Running specified simulation."
SETIDX=1
export SETIDX
sbatch -p csde -A csde slurm_batch.sbatch

