# ICER table 
# MSahu
# Update March 3, 2022

#-----------------------------------------------------------------------------------------------------------------------------

# Set up Discount Rate for Sensitivity Analysis

if (sn_dr_main == T) {DR = discount_rate[2]*100} # 3%
if (sn_dr_lo == T) {DR = discount_rate[1]*100}   # 0%
if (sn_dr_hi == T) {DR = discount_rate[3]*100}   # 5%

if (bia == T) {DR = "0"}

# Set up alternative scenario

if (alt_scen_2 == T) {SCEN = "2"}
if (alt_scen_2a == T) {SCEN = "2a"}

# Set number to round to

ROUND_cases = 0
ROUND_deaths = 0
ROUND_icers = 1
ROUND_millions = 1
ROUND_percents = 1

# -----------------------------------------------------------------------------------------------------------------------------

## GET OUTCOMES : (All Parameter Sets) - For all years up till time horizon

# Cases Averted

cases_averted <- outcomes_avrt_cum_list[[paste0("Scenario",SCEN,".HIV_incidence.dr", DR, ".main")]] %>% 
  filter(year==hrzn) %>% 
  select(-c(year, mean, min, max)) %>% 
  t() %>% as.vector()

# Deaths Averted 

deaths_averted <- outcomes_avrt_cum_list[[paste0("Scenario",SCEN,".HIV_mortality.dr", DR, ".main")]] %>% 
  filter(year==hrzn) %>% 
  select(-c(year, mean, min, max)) %>% 
  t() %>% as.vector()

# Get DALYs, QALYs, YLLs, YLDs

scen1_outcomes <- daly_list[[paste0("Scenario1.dr", DR, ".main")]] %>% 
  group_by(param_set) %>%
  filter(year <= hrzn) %>%  # time horizon
  summarize(dalys = sum(dalys),
            qalys = sum(qalys),
            ylls = sum(ylls),
            ylds = sum(ylds)) %>% 
  select(-param_set)

scen2_outcomes <- daly_list[[paste0("Scenario",SCEN,".dr", DR, ".main")]] %>% 
  group_by(param_set) %>%
  filter(year <= hrzn) %>%  # time horizon  
  summarize(dalys = sum(dalys),
            qalys = sum(qalys),
            ylls = sum(ylls),
            ylds = sum(ylds)) %>% 
  select(-param_set)

# DALYs Averted / QALYs gained
  
dalys_averted = -(scen2_outcomes$dalys - scen1_outcomes$dalys)
qalys_gained = scen2_outcomes$qalys - scen1_outcomes$qalys

# PRINT

if (ICER_TABLE == "ON") {
  
  if (DR == 0) {
  
  print("HEALTH GAINS")
    
  print(paste0("Cases Averted = ", round(mean(cases_averted), ROUND_cases), 
               " (UI = ",round(min(cases_averted), ROUND_cases) ,", ",round(max(cases_averted), ROUND_cases) ,")"))
  
  print(paste0("Deaths Averted = ", round(mean(deaths_averted), ROUND_deaths), 
               " (UI = ",round(min(deaths_averted), ROUND_deaths) ,", ",round(max(deaths_averted), ROUND_deaths) ,")"))
  
  print(paste0("DALYs Averted = ", round(mean(dalys_averted)/1e6, ROUND_millions), 
               " (UI = ",round(min(dalys_averted)/1e6, ROUND_millions) ,", ",round(max(dalys_averted)/1e6, ROUND_millions) ," million)"))
  }
  
} 

if (ICER_TABLE == "OFF") {
  
  print("HEALTH GAINS")
  
  print(paste0("Cases Averted = ", round(mean(cases_averted), ROUND_cases), 
               " (UI = ",round(min(cases_averted), ROUND_cases) ,", ",round(max(cases_averted), ROUND_cases) ,")"))
  
  print(paste0("Deaths Averted = ", round(mean(deaths_averted), ROUND_deaths), 
               " (UI = ",round(min(deaths_averted), ROUND_deaths) ,", ",round(max(deaths_averted), ROUND_deaths) ,")"))
  
  print(paste0("DALYs Averted = ", round(mean(dalys_averted)/1e6, ROUND_millions), 
               " (UI = ",round(min(dalys_averted)/1e6, ROUND_millions) ,", ",round(max(dalys_averted)/1e6, ROUND_millions) ," million)"))
}
  
# --------------------------------------------------------------------------------------------------------------------------

# GET COSTS: 

# Incremental Costs (All Parameter Sets)

inc_costs <- costs_cum_incr_list[[paste0("Scenario",SCEN,".total.dr", DR,".main")]] %>% 
  filter(year==hrzn) %>% 
  select(-c(year, mean, min, max)) %>% 
  t() %>% as.vector()

# --------------------------------------------------------------------------------------------------------------------------

# ICERs 

icer_cases = inc_costs / cases_averted   # Cost per case averted

icer_deaths = inc_costs / deaths_averted # Cost per death averted

icer = inc_costs / dalys_averted # Cost per DALY averted

icer_qaly = inc_costs / qalys_gained  # Cost per QALY gained 

# Print Discounted version for ICER table

if (ICER_TABLE == "ON") {
  
  if (DR == 3) {
    
    print("COST-EFFECTIVENESS")
    
    print(paste0("Cost per case averted, 2020 USD = ", round(mean(icer_cases), ROUND_icers), 
                 " (UI = ",round(min(icer_cases), ROUND_icers) ,", ",round(max(icer_cases), ROUND_icers) ,")"))
    
    print(paste0("Cost per death averted, 2020 USD = ", round(mean(icer_deaths), ROUND_icers), 
                 " (UI = ",round(min(icer_deaths), ROUND_icers) ,", ",round(max(icer_deaths), ROUND_icers) ,")"))
    
    print(paste0("Cost per DALY averted, 2020 USD = ", round(mean(icer), ROUND_icers), 
                 " (UI = ",round(min(icer), ROUND_icers) ,", ",round(max(icer), ROUND_icers) ,")"))
  }
    
}

if (ICER_TABLE == "OFF") {
  
  print("COST-EFFECTIVENESS")
  
  print(paste0("Cost per case averted, 2020 USD = ", round(mean(icer_cases), ROUND_icers), 
               " (UI = ",round(min(icer_cases), ROUND_icers) ,", ",round(max(icer_cases), ROUND_icers) ,")"))
  
  print(paste0("Cost per death averted, 2020 USD = ", round(mean(icer_deaths), ROUND_icers), 
               " (UI = ",round(min(icer_deaths), ROUND_icers) ,", ",round(max(icer_deaths), ROUND_icers) ,")"))
  
  print(paste0("Cost per DALY averted, 2020 USD = ", round(mean(icer), ROUND_icers), 
               " (UI = ",round(min(icer), ROUND_icers) ,", ",round(max(icer), ROUND_icers) ,")"))
}

#=============================================================================================================================

# BUDGET IMPACT ANALYSIS 

# Note: Budget Impact is always Zero Discounting

#==============================================================================================================================
 
# Import DoH Budget and inflate -----------------------------------------------------------------------------------------------

setwd("/")

library(pacman)

source('H:/repos/fgh/FUNCTIONS/currency_conversion.R')

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")

bia_params <- read.csv(paste0(cea_path,"parameters/doh_budget.csv"))  %>% 
  mutate(currency = "lcu")

doh_inflated <- currency_conversion(data = bia_params,
                                    col.value = "budget",
                                    col.loc = "iso3", 
                                    col.currency = "currency",
                                    base.year = 2020, 
                                    col.currency.year="year",
                                    base.unit = "usd",
                                    converter.version = 5.2)

doh_annual_budget = doh_inflated[doh_inflated$activity == "DOH HIV expenditure", "budget"] %>% as.numeric() # 313.2 million

# Calculate 5-year and Annual Budget Impact ------------------------------------------------------------------------------------  

bia = T  # Sets Discount Rate = 0

# Incremental Costs (All Parameter Sets)

inc_costs <- costs_cum_incr_list[[paste0("Scenario", SCEN,".total.dr", DR,".main")]] %>% 
  filter(year==hrzn) %>% 
  select(-c(year, mean, min, max)) %>% 
  t() %>% as.vector()

# Budget Impact

incr_budget <- costs_cum_incr_list[[paste0("Scenario", SCEN, ".total.dr", DR, ".main")]] %>% 
  filter(year == bia_hrzn) %>% 
  select(-c(year, mean, min, max)) %>% 
  t() %>% as.vector()

annual_budget <- incr_budget /5 

budget_impact = annual_budget / doh_annual_budget


# Annualize the cost

time_span = hrzn - 2020 + 1

# PRINT


if (ICER_TABLE == "ON") {
  
  if (DR == 0) {
    
    print("COST & BUDGET IMPACT")
    
    print(paste0("Annual incremental cost (cumulative), 2020 USD = ", round(mean(inc_costs)/1e6/time_span, ROUND_millions), " million ", 
                 " (UI = ",round(min(inc_costs)/1e6/time_span, ROUND_millions) ,", ",round(max(inc_costs)/1e6/time_span, ROUND_millions) ,")"))
    
    print(paste0("5-year Annual Program Cost, 2020 USD = ", round(mean(annual_budget)/1e6, ROUND_millions), " million ", 
                 " (UI = ",round(min(annual_budget)/1e6, ROUND_millions) ,", ",round(max(annual_budget)/1e6, ROUND_millions) ,")"))
    
    print(paste0("Community ART requires an ", round( mean(budget_impact) *100, ROUND_percents), "% initial investment",
                 " (UI = ",round( min(budget_impact) *100, ROUND_percents) ,", ",round( max(budget_impact) *100, ROUND_percents) ,")"))
  }
  
} 

if (ICER_TABLE == "OFF") {
 
  print("COST & BUDGET IMPACT")
  
  print(paste0("Total annual incremental cost, 2020 USD = ", round(mean(inc_costs)/1e6/time_span, ROUND_millions), " million ", 
               " (UI = ",round(min(inc_costs)/1e6/time_span, ROUND_millions) ,", ",round(max(inc_costs)/1e6/time_span, ROUND_millions) ,")"))
  
  print(paste0("5-year Annual Program Cost, 2020 USD = ", round(mean(annual_budget)/1e6, ROUND_millions), " million ", 
               " (UI = ",round(min(annual_budget)/1e6, ROUND_millions) ,", ",round(max(annual_budget)/1e6, ROUND_millions) ,")"))
  
  print(paste0("Community ART requires an ", round( mean(budget_impact) *100, ROUND_percents), "% initial investment",
               " (UI = ",round( min(budget_impact) *100, ROUND_percents) ,", ",round( max(budget_impact) *100, ROUND_percents) ,")"))
}

# Clean up

rm(doh_inflated, bia_params, budget_impact)
bia = F

# ---------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")