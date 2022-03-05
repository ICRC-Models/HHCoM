# MASTER SCRIPT
# Update: March 2, 2022
# MSahu

# TO DO LIST
# DALYs
# Further sensitivity analysis - assumptions

# New file for parameters "01 _ setup"

# Set up  -------------------------------------------------------------------------------------

rm(list=ls())

setwd("/")

library("tidyverse")
library("dplyr")

dir <- "C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/"
cea_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/CEA/"
main_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/DoArtOutputs/"
vmmc_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/DoArtOutputs/Sensitivity_VMMCscaleUp/"

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

hrzn <- 2060
bia_hrzn = 2024

# Threshold

threshold <- 750

# Run code ---------------------------------------------------------------------------------------

source(paste0(cea_path, "helper.R"))
source(paste0(cea_path, "01_setup_sensitivity.R"))
source(paste0(cea_path, "02_cases_deaths.R"))
source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "04_daly_qaly.R"))
source(paste0(cea_path, "05_icer_bia.R"))

# ---------------------------------------------------------------------------------------------------

