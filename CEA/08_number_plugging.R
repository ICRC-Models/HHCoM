# Number plugging script
# September 30, 2021

# --------------------------------------------------------------------------------------------------

#"and XX DALYs among men (XX%, UR: XX-XX%) 

source(paste0(cea_path, "04_daly_qaly_men.R"))

x1 <- men_daly_list[["Scenario1.dr0.main"]] %>% 
  group_by(param_set) %>%
  filter(year <= 2060) %>% # time horizon  
  summarize(dalys = sum(dalys)) 

x2 <- men_daly_list[["Scenario2.dr0.main"]]  %>% 
  group_by(param_set) %>%
  filter(year <= 2060) %>% # time horizon  
  summarize(dalys = sum(dalys)) 

diff = x2$dalys - x1$dalys
diff_mean = round(mean(diff)/1e6, 2)
diff_min = round(min(diff)/1e6, 2)
diff_max = round(max(diff)/1e6, 2)

c(diff_mean, diff_min, diff_max)

pct_diff = (x2$dalys - x1$dalys)/ x1$dalys * 100
pct_mean = round(mean(pct_diff), 1)
pct_min = round(min(pct_diff), 1)
pct_max = round(max(pct_diff), 1)

c(pct_mean, pct_min, pct_max)

# XX DALYs among women (XX%, UR: XX-XX%)."

source(paste0(cea_path, "04_daly_qaly_women.R"))

x1 <- women_daly_list[["Scenario1.dr0.main"]] %>% 
  group_by(param_set) %>%
  filter(year <= 2060) %>% # time horizon  
  summarize(dalys = sum(dalys)) 

x2 <- women_daly_list[["Scenario2.dr0.main"]]  %>% 
  group_by(param_set) %>%
  filter(year <= 2060) %>% # time horizon  
  summarize(dalys = sum(dalys)) 

diff = x2$dalys - x1$dalys
diff_mean = round(mean(diff)/1e6, 2)
diff_min = round(min(diff)/1e6, 2)
diff_max = round(max(diff)/1e6, 2)

c(diff_mean, diff_min, diff_max)

pct_diff = (x2$dalys - x1$dalys)/ x1$dalys * 100
pct_mean = round(mean(pct_diff), 1)
pct_min = round(min(pct_diff), 1)
pct_max = round(max(pct_diff), 1)

c(pct_mean, pct_min, pct_max)

# -------------------------------------------------------------------------------------------------------

# " For both Standard of Care and Home Testing + Community ART, program costs for ART delivery were projected to decline 
# between 2020 and 2060 due to reduced prevalence, however the 
# estimated decline was greater for Home Testing + Community ART (XX%, UR: XX-XX%) compared with Standard of Care (XX%, UR: XX-XX%) "

x1 <- costs_total_list[["Scenario1.total.dr0.main"]]
x2 <- costs_total_list[["Scenario2.total.dr0.main"]]

scen1_2020 = x1[x1$year == 2020, 5:29]
scen1_2060 = x1[x1$year == 2060, 5:29]
scen2_2020 = x2[x2$year == 2020, 5:29]
scen2_2060 = x2[x2$year == 2060, 5:29]

# Calculate percentage  - SCENARIO 1

pct_scen1 <- ( scen1_2060 - scen1_2020) / scen1_2020 * 100 

round(apply(pct_scen1, 1, FUN = mean), 1)
round(apply(pct_scen1, 1, FUN = min), 1)
round(apply(pct_scen1, 1, FUN = max), 1)

# Calculate percentage  - SCENARIO 2

pct_scen2 <- ( scen2_2060 - scen2_2020) / scen2_2020 * 100 

round(apply(pct_scen2, 1, FUN = mean), 1)
round(apply(pct_scen2, 1, FUN = min), 1)
round(apply(pct_scen2, 1, FUN = max), 1)

# By XX, Home Testing + Community ART was expected to be cost-saving compared with the Standard of Care"

x1 = costs_total_list[["Scenario2.total.dr0.main"]] - costs_total_list[["Scenario1.total.dr0.main"]]
x1 <- x1 %>% addYearCol() %>% recalcFuns()
year_cost_saving = min(x1[x1$mean<0, "year"])


incCosts = x1[x1$year == year_cost_saving, c("mean", "min", "max")]
incCosts 


# XX proportion fo the runs

incCosts = x1[x1$year == year_cost_saving, 5:29]

sum(incCosts <0) / 25

# -------------------------------------------------------------------------------------------------------

# "Using a shorter time horizon, 
# the estimated cost per DALY averted was XX by 2045 (UR: XX-XX) and XX by 2030 (UR: XX-XX). "

hrzn <- 2045
source(paste0(cea_path, "04_daly_qaly.R"))

-(scen2_dalys.dr0 - scen1_dalys.dr0)/1e6


