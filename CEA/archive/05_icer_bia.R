# ICER table 
# Msahu
# July 14, 2021

# TO DO
# ICERs should be calculated for each of the 25 parameter sets, and then mean min and max are taken
# CHECK with PAUL - which version of BIA to use? should the comparison be scenario 1, or DOH budget?


###########################################################################################



# DALYS 
scen1_dalys.dr0 <- NULL
scen1_dalys.dr0[1] = daly_list_summary[["Scenario1.dr0.main"]]$mean[daly_list_summary[["Scenario1.dr0.main"]]$param == "dalys"]
scen1_dalys.dr0[2] = daly_list_summary[["Scenario1.dr0.main"]]$min[daly_list_summary[["Scenario1.dr0.main"]]$param == "dalys"]
scen1_dalys.dr0[3] = daly_list_summary[["Scenario1.dr0.main"]]$max[daly_list_summary[["Scenario1.dr0.main"]]$param == "dalys"]
names(scen1_dalys.dr0) <- c("mean", "min", "max")
scen1_dalys.dr0

scen2_dalys.dr0 <- NULL
scen2_dalys.dr0[1] = daly_list_summary[["Scenario2.dr0.main"]]$mean[daly_list_summary[["Scenario2.dr0.main"]]$param == "dalys"]
scen2_dalys.dr0[2] = daly_list_summary[["Scenario2.dr0.main"]]$min[daly_list_summary[["Scenario2.dr0.main"]]$param == "dalys"]
scen2_dalys.dr0[3] = daly_list_summary[["Scenario2.dr0.main"]]$max[daly_list_summary[["Scenario2.dr0.main"]]$param == "dalys"]
names(scen2_dalys.dr0) <- c("mean", "min", "max")
scen2_dalys.dr0

-(scen2_dalys.dr0 - scen1_dalys.dr0)/1e6

(scen2_dalys.dr0 - scen1_dalys.dr0)/ scen2_dalys.dr0 * 100

scen1_dalys.dr3 <- NULL
scen1_dalys.dr3[1] = daly_list_summary[["Scenario1.dr3.main"]]$mean[daly_list_summary[["Scenario1.dr3.main"]]$param == "dalys"]
scen1_dalys.dr3[2] = daly_list_summary[["Scenario1.dr3.main"]]$min[daly_list_summary[["Scenario1.dr3.main"]]$param == "dalys"]
scen1_dalys.dr3[3] = daly_list_summary[["Scenario1.dr3.main"]]$max[daly_list_summary[["Scenario1.dr3.main"]]$param == "dalys"]
names(scen1_dalys.dr3) <- c("mean", "min", "max")
scen1_dalys.dr3

scen2_dalys.dr3 <- NULL
scen2_dalys.dr3[1] = daly_list_summary[["Scenario2.dr3.main"]]$mean[daly_list_summary[["Scenario2.dr3.main"]]$param == "dalys"]
scen2_dalys.dr3[2] = daly_list_summary[["Scenario2.dr3.main"]]$min[daly_list_summary[["Scenario2.dr3.main"]]$param == "dalys"]
scen2_dalys.dr3[3] = daly_list_summary[["Scenario2.dr3.main"]]$max[daly_list_summary[["Scenario2.dr3.main"]]$param == "dalys"]
names(scen2_dalys.dr3) <- c("mean", "min", "max")
scen2_dalys.dr3

dalys_avrt.dr3 = -(scen2_dalys.dr3 - scen1_dalys.dr3)

# QALYs

scen1_qalys.dr0 <- NULL
scen1_qalys.dr0[1] = daly_list_summary[["Scenario1.dr0.main"]]$mean[daly_list_summary[["Scenario1.dr0.main"]]$param == "qalys"]
scen1_qalys.dr0[2] = daly_list_summary[["Scenario1.dr0.main"]]$min[daly_list_summary[["Scenario1.dr0.main"]]$param == "qalys"]
scen1_qalys.dr0[3] = daly_list_summary[["Scenario1.dr0.main"]]$max[daly_list_summary[["Scenario1.dr0.main"]]$param == "qalys"]
names(scen1_qalys.dr0) <- c("mean", "min", "max")
scen1_qalys.dr0/1e6

scen2_qalys.dr0 <- NULL
scen2_qalys.dr0[1] = daly_list_summary[["Scenario2.dr0.main"]]$mean[daly_list_summary[["Scenario2.dr0.main"]]$param == "qalys"]
scen2_qalys.dr0[2] = daly_list_summary[["Scenario2.dr0.main"]]$min[daly_list_summary[["Scenario2.dr0.main"]]$param == "qalys"]
scen2_qalys.dr0[3] = daly_list_summary[["Scenario2.dr0.main"]]$max[daly_list_summary[["Scenario2.dr0.main"]]$param == "qalys"]
names(scen2_qalys.dr0) <- c("mean", "min", "max")
scen2_qalys.dr0/1e6

(scen2_qalys.dr0 - scen1_qalys.dr0)/1e6

# ------------------------------------------------------------------------------------------

# COSTS

inc_costs <- costs_cum_incr_list$Scenario2.total.dr3.main  %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

print(paste0("Total incremental cost, 2020 USD = ", round(inc_costs[1]/1e6, 2), " million ", 
             " (UI = ",round(inc_costs[2]/1e6, 2) ,", ",round(inc_costs[3]/1e6, 2) ,")"))


print(paste0("Total annual incremental cost, 2020 USD = ", round(inc_costs[1]/1e6/41, 2), " million ", 
             " (UI = ",round(inc_costs[2]/1e6/41, 2) ,", ",round(inc_costs[3]/1e6/41, 2) ,")"))


# OUTCOMES

cases_averted <- outcomes_avrt_cum_list$Scenario2.HIV_incidence.dr0.main %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

deaths_averted <- outcomes_avrt_cum_list$Scenario2.HIV_mortality.dr0.main %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

# Find total cost - percentage increase

scen1_cost <- costs_total_list$Scenario1.total.dr0.main$mean[[41]]
scen2_cost <- costs_total_list$Scenario2.total.dr0.main$mean[[41]]

pct_inc_cost = ((scen2_cost - scen1_cost)/scen1_cost)*100

#-----------------------------------------------------------------------------

# BIA

# Import DoH Budget and inflate -----------------------------------------------

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

# ----------------------------------------------------------------------------------------  

incr_budget <- costs_cum_incr_list$Scenario2.total.dr0.main %>% 
  filter(year == bia_hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix() 

annual_budget <- incr_budget /5 

print(paste0("5-year Annual Program Cost, 2020 USD = ", round(annual_budget[1]/1e6, 2), " million ", 
             " (UI = ",round(annual_budget[2]/1e6, 2) ,", ",round(annual_budget[3]/1e6, 2) ,")"))

budget_impact = annual_budget / doh_annual_budget

print(paste0("Community ART requires an ", round( budget_impact[1]*100, 1), "% initial investment",
             " (UI = ",round( budget_impact[2]*100, 1) ,", ",round( budget_impact[3]*100, 1) ,")"))

rm(doh_inflated, bia_params, budget_impact)




# ---------------------------------------------------------------------------- 

# ICER

icer <- inc_costs / dalys_avrt.dr3

print(paste0("Cost per DALY averted, 2020 USD = ", round(icer[1], 1), 
             " (UI = ",round(icer[2], 1) ,", ",round(icer[3], 1) ,")"))


# Cost per case averted

cost_case_averted <- costs_cum_incr_list$Scenario2.total.dr3.main / outcomes_avrt_cum_list$Scenario2.HIV_incidence.dr3.main
cost_case_averted$year <- 2020:2060
cost_case_averted <- cost_case_averted %>% recalcFuns()

# Cost per death averted

cost_death_averted <- costs_cum_incr_list$Scenario2.total.dr3.main / outcomes_avrt_cum_list$Scenario2.HIV_mortality.dr3.main
cost_death_averted$year <- 2020:2060
cost_death_averted <- cost_death_averted %>% recalcFuns()


# Cost per DALY averted

inc_daly = -(daly_list_summary$Scenario2.dr3.main$mean[3] - daly_list_summary$Scenario1.dr3.main$mean[3])
inc_cost = costs_cum_incr_list$Scenario2.total.dr3.main[costs_cum_incr_list$Scenario2.total.dr3.main$year==2060, 2:4]
inc_cost / inc_daly

inc_qaly = (scen2_qalys.dr0 - scen1_qalys.dr0)
inc_cost / inc_qaly


# -----------------------------------------------------------------------------

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")