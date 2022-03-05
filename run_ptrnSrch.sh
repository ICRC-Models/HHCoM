USRNAME=carajb    # SET ME: your username
echo "${USRNAME}"
export USRNAME

DIRPATH=/gscratch/csde/${USRNAME}/HHCoM    # SET ME: path to your HHCoM directory
echo "${DIRPATH}"
export DIRPATH

TCURR=0    # t_curr, iteration of calibration
echo "${TCURR}"
export TCURR

DATE=22Aug19    # SET ME: date identifier of claibration
echo "${DATE}"
export DATE

echo "Running simulations."
INT=0
for SETIDX in $(seq 1 1 2); do 
export SETIDX
sbatch -p csde -A csde --mail-user=${USRNAME}@uw.edu slurm_ptrnSrch.sbatch --qos=MaxJobs4
#INT=$(($INT + 1))
#if [ $INT -ge 5 ]; then 
#sleep 4800 # pause to be kind to the scheduler
#INT=0
#fi 
done
