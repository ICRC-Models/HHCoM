# Extract outcomes - POSTER
# July 1, 2021

hrzn <- 2060
bi_hrzn = 2024

inc_costs <- costs_cum_incr_list$Scenario2.total.dr3.main  %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

# OUTCOMES

cases_averted <- outcomes_avrt_cum_list$Scenario2.HIV_incidence.dr0.main %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

deaths_averted <- outcomes_avrt_cum_list$Scenario2.HIV_mortality.dr0.main %>% 
  filter(year==hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix()

# BIA

budget_impact <- costs_cum_incr_list$Scenario2.total.dr0.main %>% 
  filter(year == bi_hrzn) %>% 
  select(mean, min, max) %>% 
  as.matrix() 

annual_budget_impact <- budget_impact /5

# Find total cost - percentage increase

scen1_cost <- costs_total_list$Scenario1.total.dr0.main$mean[[41]]
scen2_cost <- costs_total_list$Scenario2.total.dr0.main$mean[[41]]


((scen2_cost - scen1_cost)/scen1_cost)*100

