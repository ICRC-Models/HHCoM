TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=10Dec19
echo "${DATE}"
export DATE

echo "Running MATLAB script to get matrix size."
#sbatch -p csde -A csde slurm_sizeMatrix.sbatch
#sleep 300
FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})
echo "${NSETS}" 
export NSETS

echo "Running simulations, first try."
SEQ40=($(seq 1 40 3200 ))
SEQ32=($(seq 3201 32 5760))
SEQ28p=($(seq 5761 28 8000))
SEQ28s=($(seq 8001 28 ${NSETS}))    # set up for NSETS=10240
LENGTH40=${#SEQ40[@]}    # note: all three lengths should be the same
LENGTH32=${#SEQ32[@]}
LENGTH28=${#SEQ28p[@]}
echo "${LENGTH40}"
echo "${LENGTH32}"
echo "${LENGTH28}"
for i in $(seq 1 4 ${LENGTH40}); do
    for j in $(seq $((${i}-1)) 1 $((${i}+3))); do    # submit 4 simulations for each target node at once
        SETIDX=${SEQ40[$j]}
		export SETIDX
		sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs4 --ntasks-per-node=40
        SETIDX=${SEQ32[$j]}
		export SETIDX
		sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs4 --ntasks-per-node=32
        SETIDX=${SEQ28p[$j]}
		export SETIDX
		sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs4 --ntasks-per-node=28
        SETIDX=${SEQ28s[$j]}
		export SETIDX
		sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs4 --ntasks-per-node=28
    done
	sleep 19200    # give submitted simulations time to finish 
done

#echo "Running MATLAB script to identify failed simulations."
#sbatch -p ckpt -A csde-ckpt slurm_idMissing.sbatch
#sleep 300
#FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
#RERUN=$(<$FILE)
#echo "$RERUN"

#echo "Re-running failed simulations."
#while [ ! -z "$RERUN" ]; do
#MISSING=($RERUN)
#INT=0
#for SETIDX in ${MISSING[@]}; do
#export SETIDX
#sbatch -p ckpt -A csde-ckpt slurm_batch.sbatch --qos=MaxJobs4
#INT=$(($INT + 1))
#if [ $INT -ge 50 ] || [ $INT -eq ${#MISSING[@]} ]; then 
#sleep 5 #4800 # pause to be kind to the scheduler
#INT=0
#fi 
#done

#echo "Running MATLAB script to identify failed simulations, again."
#sbatch -p ckpt -A csde-ckpt slurm_idMissing.sbatch
#sleep 300
#FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
#RERUN=$(<$FILE)
#echo "$RERUN"
#done

#echo "Running MATLAB abc_smc script to get next set of particles."
#sbatch -p csde -A csde slurm_abc.sbatch
#sleep 21600
 
#echo "Running MATLAB idParamRanges script to get ranges of parameters in best-fitting sets."
#sbatch -p ckpt -A ckpt-csde slurm_idParamRanges.sbatch

