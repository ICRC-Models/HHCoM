TCURR=0    # t_curr
echo "${TCURR}"
export TCURR

DATE=19Dec19
echo "${DATE}"
export DATE

CUTOFF=6132

tail -n 10020 ./Params/orderedLL_calib_${DATE}_${TCURR}.dat |
while IFS=',' read a b; do
  echo "$a"
done > "./Params/weededLL_calib_${DATE}_${TCURR}.dat"

FILE=./Params/weededLL_calib_${DATE}_${TCURR}.dat
SETS=$(<${FILE})
#echo "${SETS}" 

for i in ${SETS}; do
    if [ ${TCURR} == 0 ]; then
	    rm ./HHCoM_Results/toNow_${DATE}_noBaseVax_baseScreen_hpvHIVcalib_${TCURR}_${i}.mat
	fi
	if [ ${TCURR} -gt 0 ]; then
	    if [ ${i} -gt ${CUTOFF} ]; then
		    NEWI=$((${i}-${CUTOFF}))
		    rm ./HHCoM_Results/toNow_${DATE}_noBaseVax_baseScreen_hpvHIVcalib_${TCURR}_${NEWI}.mat
		fi
	fi
done

