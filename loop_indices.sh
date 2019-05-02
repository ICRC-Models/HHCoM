for PARAM1 in $(seq 1 5); do 
# 
echo "${PARAM1}" 
export PARAM1
# 
sbatch -o out_p${PARAM1}.stdout.txt \ 
-e out_${PARAM1}.stdout.txt \ 
--job-name=jobArrayTest_${PARAM1} \ 
vary_params.sbatch 
# 
sleep 1 # pause to be kind to the scheduler 
done
