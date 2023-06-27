The Kenya_treatment branch uses Gui's original Kenya_National branch and adds in the following components: 

- Single dose vaccine efficacy from the KEN-SHE results presented at IPVC 2023
- CC treatment and symptomatic detection
- Gradual scale up of vaccine coverage 

This branch was developed for the Kenya single dose CEA. 

Author: Christine L Hathaway 
Date: May 24, 2023

Workflow for single dose CEA: 
- historicalSim and futureSim 
- Run calib2_sumll4sets on the cluster. Make sure to comment/uncomment historical or future sim depending on which one you are running. 
- For MATLAB to CSV, use vaxCEA_multSims_mainFunction, which calls vaxCEA_multSims_processResults 
- In KECEA, cleanMatlabOutputs.R transforms CSV files into CSV template files for the health econ analysis
- modelAnalysis.Rmd transforms CSV files to produce figures

Other notes: 
- Scenario 10 is the only one that has a more condensed version of ccSymp and ccTreat. I did this because the cluster was running out of memory. So that means scenario 10 must be processed different from the other scenarios. 