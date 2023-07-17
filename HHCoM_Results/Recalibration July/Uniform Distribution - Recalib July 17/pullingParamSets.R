# title: Pulling parameter sets 
# author: Christine L Hathaway
# date: July 17, 2023
# description: For the recalibration of kSymp, this was me pulling from uniform distributions and saving all the parameter sets in an excel file to be ready by sympCalibration.m

library(readxl)
library(tidyverse)
library(lubridate)
library(janitor)
library(tidyr)
library(ggplot2)
library(openxlsx)

kSymp_L <- runif(100, 0, 0.03 )
kSymp_R <- runif(100, 0, 0.8)
kSymp_D <- runif(100, 0.7, 1)

paramSet <- data.frame("L" = kSymp_L, "R" = kSymp_R, "D" = kSymp_D)

OUT <- createWorkbook()

addWorksheet(OUT, "Param sets")

writeData(OUT, x = paramSet, startRow = 1, sheet="Param sets")

saveWorkbook(OUT, "/Users/clh89/MATLAB/Projects/Kenya_treatment/Params/kSymp_ParamSets.xlsx", overwrite = TRUE)