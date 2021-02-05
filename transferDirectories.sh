cd HHCoM_Results/

simList=( 1 2 3 6 8 9 11 12 13 15 20 21 22 26 27 32 34 35 38 39 40 41 42 45 47 )
for i in ${simList[@]}; do
    mkdir Vaccine22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_discontFxd_WHO-SCES34_6_$i
	cd Vaccine22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_discontFxd_WHO-SCES34_6_$i
	scp carajb@mox.hyak.uw.edu:/gscratch/csde/carajb/HHCoM/HHCoM_Results/Vaccine22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_discontFxd_WHO-SCES34_6_$i/vaxSimResult*.mat .
	cd ..
done
