NSETS=3088    # nSets
echo "${NSETS}" 
export NSETS

INT=0

#for SETIDX in $(seq 1 16 ${NSETS}); do 
MISSING=(529 609 625 689 785 849 945 1009 1169 1265 1297 1329 1361 1425 1489 1585 1649 1745 2273 2369 2385 2929 3009 3025 3089 3105)

for SETIDX in ${MISSING[@]}; do
export SETIDX

INT=$(($INT + 1))

sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs11

if [ $INT -ge 10 ]; then 
sleep 6000 # pause to be kind to the scheduler
INT=0
fi 

done
