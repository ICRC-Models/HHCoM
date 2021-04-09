TCURR=1    # t_curr
echo "${TCURR}"
export TCURR

DATE=26Oct20
echo "${DATE}"
export DATE

echo "Running specified simulation."
SETIDX=1
export SETIDX
sbatch -p csde -A csde slurm_runMultHistSims.sbatch

