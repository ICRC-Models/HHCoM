cd HHCoM_Results/

for i in $(seq 1 1 25); do
    mkdir Vaccine22Apr20Ph2V2_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd_diagHiv075_DoART_S2_11_$i
	cd Vaccine22Apr20Ph2V2_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd_diagHiv075_DoART_S2_11_$i
	scp carajb@mox.hyak.uw.edu:/gscratch/csde/carajb/HHCoM/HHCoM_Results/Vaccine22Apr20Ph2V2_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd_diagHiv075_DoART_S2_11_$i/*.mat .
	cd ..
done
