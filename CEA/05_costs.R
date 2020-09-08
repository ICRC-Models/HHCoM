# Generate COST over time per capita for each Scenario
# MSahu
# September 8, 2020

library(dplyr)
library(pacman)
source("00_master.R")
source('H:/repos/fgh/FGH_2019/04_functions/currency_conversion.R')

# TO DO:
# 1. Inflate any costs -- DONE
# 2. Discount
# 3. Double check if reasonable to assume everyone on ART is CD4>350
# 4. Check ART costs
# 5. include min/max values as well
# 6. Include startup cost cb_art , 
# 7. Cumulative costs not needed
# 8. Calculate QALYS!!

# Import DF with literature costs and inflate

cost_params <- read.csv(paste0(cea_path,"costs/cost_params.csv")) %>% 
  mutate(iso3="USA") 

cost_inflated <- currency_conversion(data = cost_params,
                      col.value = "param_val",
                      col.loc = "iso3", 
                      base.year = 2020, 
                      col.currency="currency", 
                      col.currency.year="year")

# Import prevalence DFs




# Annual costs, 2030 - Scenario 1

ART_cost <- prev_ART %>% 
  select(year, mean1) %>% 
  mutate(ART_cost = mean1*art + mean1*hosp_art, # ART cost plus hospitalization cost
         cum_ART_cost = cumsum(ART_cost))

cost_200 <- prev_200 %>% 
  select(year, mean1) %>% 
  mutate(hosp200_cost = mean1*hosp_200,
         cum_hosp200_cost = cumsum(hosp200_cost))

cost_200_350 <- prev_200_350 %>% 
  select(year, mean1) %>% 
  mutate(hosp200_350_cost = mean1*hosp_200_350,
         cum_hosp200_350_cost = cumsum(hosp200_350_cost))

cost_350  <- prev_350 %>% 
  select(year, mean1) %>% 
  mutate(hosp350_cost = mean1*hosp_350,
         cum_hosp350_cost = cumsum(hosp350_cost))

total_cost1 <- ART_cost %>% 
  left_join(cost_200, by="year") %>% 
  left_join(cost_200_350, by="year") %>% 
  left_join(cost_350, by="year") %>% 
  mutate(total_cost1 = cum_ART_cost+cum_hosp200_cost+cum_hosp200_350_cost+cum_hosp350_cost) %>% 
  select(year, total_cost1)

# Scenario 2 cost per person  -- includes testing costs

cost_testing <- data.frame("year"=2017:2060) %>% 
  mutate(test_cost=ifelse(year%%5==0,testing_pc,0)) %>% 
  mutate(cost_test_cum=cumsum(test_cost))

ART_cost <- prev_ART %>% 
  select(year, mean2) %>% 
  mutate(ART_cost = mean2*art*(51/70) + mean2*cb_art*(19/70) + mean2*hosp_art, # ART cost plus hospitalization cost
         cum_ART_cost = cumsum(ART_cost))

cost_200 <- prev_200 %>% 
  select(year, mean2) %>% 
  mutate(hosp200_cost = mean2*hosp_200,
         cum_hosp200_cost = cumsum(hosp200_cost))

cost_200_350 <- prev_200_350 %>% 
  select(year, mean2) %>% 
  mutate(hosp200_350_cost = mean2*hosp_200_350,
         cum_hosp200_350_cost = cumsum(hosp200_350_cost))

cost_350  <- prev_350 %>% 
  select(year, mean2) %>% 
  mutate(hosp350_cost = mean2*hosp_350,
         cum_hosp350_cost = cumsum(hosp350_cost))

total_cost2 <- ART_cost %>% 
  left_join(cost_200, by="year") %>% 
  left_join(cost_200_350, by="year") %>% 
  left_join(cost_350, by="year") %>% 
  left_join(cost_testing, by="year") %>% 
  mutate(cost_test_cum=ifelse(is.na(cost_test_cum),0,cost_test_cum)) %>%  #replace NA with zero
  mutate(total_cost2 = cum_ART_cost+cum_hosp200_cost+cum_hosp200_350_cost+cum_hosp350_cost+cost_test_cum) %>% 
  select(year, total_cost2)


cost_diff <- total_cost1 %>% 
  left_join(total_cost2, by="year") %>% 
  mutate(cost_diff_1_2=total_cost2-total_cost1)

