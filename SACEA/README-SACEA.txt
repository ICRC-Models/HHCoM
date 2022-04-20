README-SACEA.txt

LAST UPDATED: 4/15/2022
WRITTEN BY: Christine Hathaway 

DESCRIPTION: For the CEA of the South Africa screening paper, Matlab code was written to pull results from the model for the 10 scenarios (S0-S9), over 25 parameters. The following files are required for pulling these results: 

In the main directory, `vaxCEA_multSims_SACEA_CH.m`
In the main directory, `loopingCEAOverScenarios_v2.m`
In the SACEA folder, `cleanMatlabOutputs.R`
In the main directory, `vaxCEA_multSims_CIs_SACEA.m` -- this is the code Cara wrote for the screening/treatment portion. 

To run the code, open `loopingCEAOverScenarios_v2.m`. 

Things to check:
- In `loopingCEAOverSceanrios_v2.m`, check that the fileInds are the same parameters you want to use for the CEA 
- In `vaxCEA_multSims_SACEA_CH.m`, check lastYear is the correct last year you want to run for the simulation

Go ahead and run `loopingCEAOverSceanrios_v2.m`. It loops through all the scenarios and calls the vaxCEA_multSims_SACEA function to add results for each parameter into a master result matrix. Because each loop is dependent on the last, I did not parallelize it so it may be a little slow. You can keep track of progress by a display that pops up after each scenario/parameter combination has finished running. 

It saves 6 files in the SACEA directory: 
deaths.csv
screenTreat.csv
hpvHealthSate.csv
hivHealthState.csv
totalPerAge.csv
vax.csv

Next, within the SACEA directory, open `clanMatlabOutputs.R` in RStudio. Note that I call the tidyverse, janitor, ggplot2, and lubridate libraries, but I only really use tidyverse in this code. 

Double check that paramCateg is the correct parameters list that matches up with the Matlab code. 

Adjust the working directory. 

This code reads in all the CSV files and cleans them in the same format Jacinda requested for the CEA. Each file is outputted to SACEA/Outputs, with one csv file per scenario. 

Some notes about the results: 
allDeath, nonavalCount, and bivalCount are not stratified by age, so results can be found under the age "all.ages". 