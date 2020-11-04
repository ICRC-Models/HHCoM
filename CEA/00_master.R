# MASTER SCRIPT
# September 7, 2020
#MSahu

# Set up 

rm(list=ls())

library("dplyr")

main_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/DoArtOutputs/"
cea_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/CEA/"

# Set up helper discount functions

discount_rate <- .03   
discount <- function(year_discount,discount_rate) {
  1/((1 + discount_rate)^year_discount)
}
discounter <- function(x,discount_amt) {
  x*discount_amt
}


# Set up DO ART parameters: enrollment assumptions

females_ART_scen1 <- .6223
females_ART_scen2 <- .657
females_ART_scen3 <- .857

males_ART_scen1 <- .3978
males_ART_scen2 <- .655
males_ART_scen3 <- .857

males_enrolment_clinic <- 51/72
males_enrolment_cbART <- 21/72

females_enrolment_clinic <- 70/73
females_enrolment_cbART <- 3/73

# Set up DO ART % tested scalar

DOARTpct_tested <- 0.9

# SENSITIVITY ANALYSIS  - VMMC SCALE-UP SCENARIO

main_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/DoArtOutputs/Sensitivity_vmmcScaleUpSA/"

source("01_cases_averted.R")
source("02_deaths_averted.R")
source("03_QALYS_gained.R")
source("04_DALYs_averted.R")
source("05_costs.R")