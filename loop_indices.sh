TMIN=0
TMAX=0


for ((iter = $TMIN; iter <= $TMAX; iter++))
do
    TCURR=$iter    # calibration iteration
    echo "${TCURR}"
    export TCURR

    DATE=28Jun21    # date
    echo "${DATE}"
    export DATE

    # Notes on workflow: comment out sections for previous steps when you move on to the next step
    # using    : <<'END'    at the beginning of the section and    END    at the end of the section

    # STEP ONE - calculate and save nsets in a file
    # Notes: slurm_sizeMatrix.sbatch must successfully run and create matrixSize_calib_[date]_[t_curr].dat 
    # before you can successfully run slurm_batch.sbatch

    echo "Running MATLAB script to get matrix size."
    sbatch -p csde -A csde slurm_sizeMatrix.sbatch
    sleep 180



    # STEP TWO - run simulations with potential parameters sets in parallel
    FILE=./Params/matrixSize_calib_${DATE}_${TCURR}.dat
    NSETS=$(<${FILE})
    echo "${NSETS}" 
    export NSETS
    echo "Running simulations, first try."
    SEQ28all=($(seq 1 28 ${NSETS}))    # set up for NSETS=5600
    LENGTH28=${#SEQ28all[@]}
    echo "${LENGTH28}"
    NSLKEY=$(($NSETS * 9 / 7))


    for i in $(seq 1 5 ${LENGTH28}); do
        for j in $(seq $((${i}-1)) 1 $((${i}+3))); do    # submit batches of 4 jobs (28 simulations per job) to not overwhelm scheduler
            SETIDX=${SEQ28all[$j]}
                export SETIDX
                sbatch -p ckpt -A csde-ckpt slurm_batch.sbatch
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

    HOURS=60 ##max number of hours to check before giving up

    while [ $HOURS -gt 0 ]
    do
        
        jobQ=$(squeue -u willmin)


        echo "JOBQ length ${#jobQ}"

        if [ ${#jobQ} -eq 84 ]; then
            if [ -f ./Params/negSumLogL_calib_${DATE}_${TCURR}.dat ]; then
                NegsumlogLines=$( wc -l < ./Params/negSumLogL_calib_${DATE}_${TCURR}.dat )
                if [[ $NegsumlogLines -eq $NSLKEY ]]; then
                    echo "Running MATLAB abc_smc script to get next set of particles."
                    sbatch -p ckpt -A csde-ckpt slurm_abc.sbatch
                    HOURS=0

                else
                    HOURS=60
                    # STEP THREE - check for failed simulations
                    echo "Running MATLAB script to identify failed simulations."
                    sbatch -p csde -A csde slurm_idMissing.sbatch
                    sleep 300


                    # STEP FOUR - rerun failed simulations

                    FILE=./Params/missingSets_calib_${DATE}_${TCURR}.dat
                    RERUN=$(<$FILE)
                    echo "$RERUN"
                    echo "Re-running failed simulations."
                    #while [ ! -z "$RERUN" ]; do
                        MISSING=($RERUN)
                        INT=0
                        for i in $(seq 1 5 ${LENGTH28}); do
                            for j in $(seq $((${i}-1)) 1 $((${i}+3))); do    # submit batches of 4 jobs (28 simulations per job) to not overwhelm scheduler
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
                fi 

            else
                echo "error, negsumlog file is not readable"
                HOURS=0
            fi
        fi
        HOURS=$(($HOURS - 1))
        sleep 3600
    done
done


