# MASTER SCRIPT
# September 7, 2020
# MSahu

# Set up 

rm(list=ls())

library("dplyr")

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")
main_path <- "DoArtOutputs/"
cea_path <- "CEA/"
helper_path <- "CEA/"

# PRIMARY ANALYSIS

source(paste0(helper_path, "helper.R"))
source(paste0(cea_path, "01_cases_averted.R"))
source(paste0(cea_path, "02_deaths_averted.R"))
source(paste0(cea_path, "03_QALYS_gained.R"))
source(paste0(cea_path, "04_DALYs_averted_dr0.R"))
source(paste0(cea_path, "05_costs_dr3.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "07_ICER_table.R"))  
source(paste0(helper_path, "helper.R"))


# SENSITIVITY ANALYSIS  - VMMC SCALE-UP SCENARIO

main_path <- "DoArtOutputs/Sensitivity_VMMCscaleUp/"
cea_path <- "CEA/Sensitivity_VMMCscaleUp/"

source(paste0(helper_path, "helper.R"))
source(paste0("CEA/01_cases_averted.R"))
source(paste0("CEA/02_deaths_averted.R"))
source(paste0("CEA/03_QALYS_gained.R"))
source(paste0("CEA/04_DALYs_averted.R"))
# source(paste0("CEA/05_costs_VMMC.R"))  # MUST BE CONNECTED TO VPN, or will get error

# SENSITIVITY ANALYSIS  
# 1. Local life expectancy 
#2. VARY TIME HORIZON (2030, 2045)


# To Do List
# 1. Fix DALYS - discounting and all-cause mortality. Also fix discontinuation/continuation
# 2. Fix ICER table for cleaner output
# Costs: add hospitalization costs for incident cases?
# 3. Fix NHB tables to neatly output
# 4. Allow for variation of costs - lower and upper bounds
# 5. Tornado plot!!
# 6. Scenario 3 analysis
# 7. Scenario 2a analysis/
# 8. Check undiscounted numbers match Cara's
