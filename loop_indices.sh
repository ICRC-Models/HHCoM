NSETS = 48    # nSets
echo "${NSETS}" 
export NSETS

for SETIDX in $(seq 1 16 NSETS); do 
echo "${SETIDX}" 
export SETIDX

sbatch -p csde -A csde slurm_batch.sbatch
 
sleep 1 # pause to be kind to the scheduler 

done