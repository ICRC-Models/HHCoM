TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE="05Aug19"
echo "${DATE}"
export DATE

echo "Running MATLAB script to get matrix size."
srun -p csde -A csde --time=0:10:00 --mem=58G --pty /bin/bash
export TZ="America/Los_Angeles"
module load matlab_2018a
matlab -nodisplay -nosplash -r "idSizeMatrix(${TCURR} , ${DATE})"
sleep 670
FILE="matrixSize_calib_$DATE_$TCURR.dat"
NSETS=$(<$FILE)
#NSETS=3088    # nSets
echo "${NSETS}" 
export NSETS

echo "Running simulations, first try."
INT=0
for SETIDX in $(seq 1 16 ${NSETS}); do 
export SETIDX
INT=$(($INT + 1))
sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs11
if [ $INT -ge 10 ]; then 
sleep 6000 # pause to be kind to the scheduler
INT=0
fi 
done

echo "Running MATLAB script to identify failed simulations."
srun -p csde -A csde --time=00:10:00 --mem=58G --pty /bin/bash
export TZ="America/Los_Angeles"
module load matlab_2018a
matlab -nodisplay -nosplash -r "idMissingSets(${TCURR} , ${DATE} , ${NSETS})"
sleep 670
FILE="missingSets_calib_$DATE_$TCURR.dat"
RERUN=$(<$FILE)
echo "$RERUN"

echo "Re-running failed simulations."
while [-z "$RERUN"]; do
MISSING=($RERUN)
INT=0
for SETIDX in ${MISSING[@]}; do
export SETIDX
INT=$(($INT + 1))
sbatch -p csde -A csde slurm_batch.sbatch --qos=MaxJobs11
if [ $INT -ge 10 ]; then 
sleep 6000 # pause to be kind to the scheduler
INT=0
fi 
done

echo "Running MATLAB script to identify failed simulations, again."
srun -p csde -A csde --time=00:10:00 --mem=58G --pty /bin/bash
export TZ="America/Los_Angeles"
module load matlab_2018a
matlab -nodisplay -nosplash -r "idMissingSets(${TCURR} , ${DATE} , ${NSETS})"
sleep 670
FILE="missingSets_calib_$DATE_$TCURR.dat"
RERUN=$(<$FILE)
echo "$RERUN"
done

echo "Running MATLAB abc_smc script to get next set of particles."
srun -p csde -A csde --time=04:00:00 --mem=58G --pty /bin/bash
export TZ="America/Los_Angeles"
module load matlab_2018a
matlab -nodisplay -nosplash -r "abc_smc(${TCURR} , ${DATE} , ${NSETS})"
