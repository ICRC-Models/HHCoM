TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=22Apr20Ph2V11
echo "${DATE}"
export DATE

echo "Running MATLAB script to get matrix size."
sbatch -p csde -A csde slurm_sizeMatrix.sbatch
sleep 180
FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
NSETS=$(<${FILE})
echo "${NSETS}" 
export NSETS

echo "Running simulations, first try."
SEQ28all=($(seq 1 28 ${NSETS}))      # set up for NSETS=7644
LENGTH28=${#SEQ28all[@]}
echo "${LENGTH28}"
for i in $(seq 1 5 ${LENGTH28}); do
    for j in $(seq $((${i}-1)) 1 $((${i}+3))); do    # submit 4 simulations for each target node at once 
        SETIDX=${SEQ28all[$j]}
            export SETIDX
            sbatch -p csde -A csde slurm_batch.sbatch #--qos=MaxJobs4
    done
    
    sleep 10
    if [ $i -ge 73 ] && [ $i -lt 77 ]; then
        echo "${i}"
        sleep 5400    # give submitted simulations time to finish 
    fi
    if [ $i -ge 154 ] && [ $i -lt 158 ]; then
        echo "${i}"
        sleep 5400
    fi
    if [ $i -ge 234 ] && [ $i -lt 238 ]; then
        echo "${i}"
        sleep 5400
    fi
done

: <<'END'
echo "Running MATLAB script to identify failed simulations."
sbatch -p ckpt -A csde-ckpt slurm_idMissing.sbatch
#sleep 300
#FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
#RERUN=$(<$FILE)
#echo "$RERUN"
#: <<'END'
echo "Re-running failed simulations."
#while [ ! -z "$RERUN" ]; do
    MISSING=($RERUN)
    INT=0
    for i in $(seq 1 5 ${LENGTH28}); do
	for j in $(seq $((${i}-1)) 1 $((${i}+3))); do    # submit 4 simulations for each target node at once
             SETIDX=${SEQ28all[$j]}
             if [[ " ${MISSING[@]} " =~ " ${SETIDX} " ]]; then
                 export SETIDX
                 sbatch -p ckpt -A csde-ckpt slurm_batch.sbatch #--qos=MaxJobs4
                 INT=$(($INT + 1))
             fi
        done
	if [ $INT -ge 50 ]; then 
            sleep 7200 #10800    # give submitted simulations time to finish 
	    INT=0
	fi
    done
    #if [ $INT -gt 0 ] && [ $INT -lt 5 ]; then
    #    sleep 10800    # give submitted simulations time to finish 
    #fi

    #echo "Running MATLAB script to identify failed simulations, again."
    #sbatch -p csde -A csde slurm_idMissing.sbatch
    #sleep 300
    #FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
    #RERUN=$(<$FILE)
    #echo "$RERUN"
#done
END

#echo "Running MATLAB abc_smc script to get next set of particles."
#sbatch -p ckpt -A csde-ckpt slurm_abc.sbatch
#sleep 21600
 
#echo "Running MATLAB idParamRanges script to get ranges of parameters in best-fitting sets."
#sbatch -p csde -A csde slurm_idParamRanges.sbatch


