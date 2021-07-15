# MASTER SCRIPT
# Update: July 21, 2021
# MSahu

# To Do List
# 1. Fix DALYS - discounting and all-cause mortality. Also fix discontinuation/continuation
# 2. Fix ICER table for cleaner output

# 3. Fix NHB tables to neatly output
# 5. Tornado plot!!

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

hrzn <- 2060
bia_hrzn = 2024


## SENSITIVITY ANALYSIS SETUP (turn off) -------------------------------------------------------------

# Set up DF with list of sensitivitity analyses

lengthSN = 5 # Number of sensitivity analyses

snDF <- data.frame(sn_type = rep("Cost", lengthSN),
                   sn_name = c("Home testing", "Hospitalization", "Community ART", "Clinic ART", "VMMC"),
                   sn_no = c(1:lengthSN))

# Bounds

bound <- c("lower", "upper")
snDF <- snDF %>% crossing(bound) %>% mutate(bound_abb = ifelse(bound == "lower", "lb", "ub")) %>% 
  arrange(sn_no, bound)

# Labels

snDF <- snDF %>% 
  mutate(name_abb = paste0( "sn" , sn_no, ".", bound_abb),
         sn_name_full = paste(sn_type, "for", sn_name, " - ", bound, sep = " "))

# Set on or off!

snDF$ON_or_OFF = F

# Primary settings

vs_scalar = T # include the scalar to get from ART + VS to ART
vmmc_cost = F # do not include VMMC costs

# ---------------------------------------------------------------------------------------------------

source(paste0(cea_path, "helper.R"))
source(paste0(cea_path, "01_cases_deaths.R"))
source(paste0(cea_path, "02_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "04_icer_bia.R"))

# TO DO
# source(paste0(cea_path, "03_daly_qaly.R"))
# source(paste0(cea_path, "04_ICER_table.R"))  


# ---------------------------------------------------------------------------------------------------
