TCURR=11    # t_curr
echo "${TCURR}"
export TCURR

DATE=22Apr20Ph2V2
echo "${DATE}"
export DATE

echo "Get matrix size."
FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})
echo "${NSETS}" 
export NSETS

echo "Running specified simulation."
for i in $(seq 1 5 25); do
    SETIDX=$i  #1
    export SETIDX
    sbatch -p csde -A csde slurm_runMultSims.sbatch
done

