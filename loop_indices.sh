NSETS=3110    # nSets
echo "${NSETS}" 
export NSETS

INT=0

for SETIDX in $(seq 1 16 ${NSETS}); do 
#MISSING=(2145 2289 2545 2689 3057 3105)
#for SETIDX in ${MISSING[@]}; do
export SETIDX

INT=$(($INT + 1))

sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs11

if [ $INT -ge 10 ]; then 
sleep 6000 # pause to be kind to the scheduler
INT=0
fi 

done
