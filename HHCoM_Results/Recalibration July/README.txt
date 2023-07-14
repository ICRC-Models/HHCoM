Author: Christine L Hathaway
Date created: July 13, 2023

Description: 
Everything in this folder contains the recalibration done to kSymp and the L to R, and R to D progression probabilities in the Kenya model with treatment. 

Order of running code: 
1) Uses sympCalibration (in main directory), which calls modifiedhistoricalSim (in main directory)
2) A list of tested out symp params are in this directory under sympParams.m
3) You save all the outputs from the cluster into this folder and run processRecalibResults.m, saves all data into a csv file. 
4) You can visualize/process results using vizRecalibResults.Rmd