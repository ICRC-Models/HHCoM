# MASTER SCRIPT
# Update: July 21, 2021
# MSahu

# Set up  -------------------------------------------------------------------------------------

rm(list=ls())

library("dplyr")

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")
main_path <- "DoArtOutputs/"
vmmc_path <- "DoArtOutputs/Sensitivity_VMMCscaleUp/"
cea_path <- "CEA/"

# Set up parameters ---------------------------------------------------------------------------

# Scenarios

scenarios <- c("Scenario1", "Scenario2", "Scenario2a")
names(scenarios) <- c("Standard of care",
                      "HTC + Community ART",
                      "HTC only")
alt_scenarios <- scenarios[!scenarios =="Scenario1"]

# VMMC Scale-Up

version <- c(main_path, vmmc_path)
names(version) <- c("main","vmmc")

# Discount rate of 3%

discount_rate <- c(0, .03, .05)
names(discount_rate) <- c("dr0", "dr3", "dr5")

# Column names 

df_names <- c("year","mean", "min","max", paste0("s",1:25))

# Time horizon

horizon_year <- 2060


# PRIMARY ANALYSIS -----------------------------------------------------------------------------

# Primary settings

vs_scalar = "on" # include the scalar to get from ART + VS to ART
vmmc_cost = "off" # do not include VMMC costs

source(paste0(cea_path, "helper.R"))
source(paste0(cea_path, "01_cases_deaths.R"))
source(paste0(cea_path, "05_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error

# source(paste0(cea_path, "03_QALYS_gained.R"))
source(paste0(cea_path, "04_DALYs_averted.R"))

source(paste0(cea_path, "07_ICER_table.R"))  
source(paste0(cea_path, "helper.R"))

## SENSITIVITY ANALYSIS  
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
