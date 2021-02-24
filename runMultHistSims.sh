TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=14SEP
echo "${DATE}"
export DATE

echo "Running specified simulation."
SETIDX=1
export SETIDX
sbatch -p ckpt -A csde-ckpt slurm_batch.sbatch

