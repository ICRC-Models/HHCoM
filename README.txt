The Kenya_treatment branch uses Gui's original Kenya_National branch and adds in the following components: 

- Single dose vaccine efficacy from the KEN-SHE results presented at IPVC 2023
- CC treatment and symptomatic detection
- Gradual scale up of vaccine coverage 
- Vaccination Protection of HPV + Females 

This branch was developed for the Kenya HPV Prevalanet Vaccination. 

Author: Sydney Klein 
Date: December 9th, 2024

Workflow for HPV + Vaccination: 
- historicalSim and futureSim 
- Run calib2_sumll4sets on the cluster. Make sure to comment/uncomment historical or future sim depending on which one you are running. 
- In futureSim, make sure to change loadUp to loadup2_NoIttEff or loadUp2_IttEff depending on what you are running
    - loadup2_NoIttEff has no vaccine efficacy for those with HPV infections at time of vaccination
    - loadUp2_IttEff has vaccine efficacy for those with HPV infections at time of vaccination, with data coming from the KEN SHE study
    - vaccination of hpv + individuals is accounted for in hpvVaxCU and hpvVaxSchool
- For MATLAB to CSV, use vaxCEA_multSims_mainFunction, which calls vaxCEA_multSims_processResults 


Other notes: 
- When running the matlab to csv code for the waning scenarios, make sure that waning = 0 and singleDoseBool = 0 in loadup2
- Calibration of kSymp used sympCalibration, debugTreatment_sympCalibration_historicalSim, and debugTreatment_sympCalibration_historicalSim_runthis. 

Recalibration (Jul 10, 2023, completed Jul 20, 2023): 
- Gui originally use Campos progression probabilities which are monthly. We need yearly. I adjusted prog
probs to align with Cari and Evan models. 
- Everything related to this can be found within HHCoM Results. 
- I ran sympCalibration, which calls modifiedhistoricalSim. I realized you do not have to save the full results file. We only care about the stage distribution at one point in time and the counts of CC cases. 
- Process results using processRecalibResults. 
- Then you can read in the subsequent CSV into R to filter through the results. 

Running the model with BOTH waning and catch-up: 
- I edited mixInfect in order to allow for this. Note that I coded it only specific to catch-up of age 10-19. 
- You will want to first figure out what your catch-up coverage will be. Then run the model locally and fill out 
"Adding catch-up and waning.xlsx". Input the first 6 time steps as vaxCUNormRatio. 
- Make sure waning = 1 and vaxCU = 1.