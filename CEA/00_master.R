# MASTER SCRIPT
# September 7, 2020
# MSahu

# Set up 

rm(list=ls())

library("dplyr")

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/")
main_path <- "DoArtOutputs/"
cea_path <- "CEA/"
helper_path <- "CEA/"

source(paste0(helper_path, "helper.R"))
source(paste0(cea_path, "01_cases_averted.R"))
source(paste0(cea_path, "02_deaths_averted.R"))
source(paste0(cea_path, "03_QALYS_gained.R"))
source(paste0(cea_path, "04_DALYs_averted.R"))
source(paste0(cea_path, "05_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error


# SENSITIVITY ANALYSIS  - VMMC SCALE-UP SCENARIO

main_path <- "DoArtOutputs/Sensitivity_VMMCscaleUp/"
cea_path <- "CEA/Sensitivity_VMMCscaleUp/"

source(paste0(helper_path, "helper.R"))
source(paste0("CEA/01_cases_averted.R"))
source(paste0("CEA/02_deaths_averted.R"))
source(paste0("CEA/03_QALYS_gained.R"))
source(paste0("CEA/04_DALYs_averted.R"))
source(paste0("CEA/05_costs_VMMC.R"))  # MUST BE CONNECTED TO VPN, or will get error

# SENSITIVITY ANALYSIS  
# 1. Local life expectancy 
#2. VARY TIME HORIZON (2030, 2045)


# To Do List
# 1. Fix DALYS - discounting and all-cause mortality
# 2. Fix ICER table for cleaner output
# Costs: add hospitalization costs for incident cases?
# 3. Fix NHB tables to neatly output
# 4. Allow for variation of costs - lower and upper bounds
# 5. Tornado plot!!
# 6. Scenario 3 analysis
# 7. Check undiscounted numbers match Cara's
