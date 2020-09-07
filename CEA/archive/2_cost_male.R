# Generate COST over time per capita for each Scenario
# MSahu
# July 28, 2020

library(dplyr)


# TO DO:
# 1. Inflate any costs
# 2. Discount
# 3. Double check if reasonable to assume everyone on ART is CD4>350
# 4. Check ART costs
# 5. include min/max values as well
# 6. Include startup cost cb_art

# LITERATURE COSTS

testing_pc <- 8.58 # Sharma et al, Home HTC (need to inflate from 2012 USD to 2020)

hosp_200 <- 121
hosp_200_350 <- 58
hosp_350 <- 39
hosp_art <- 45  # This assumes everyone on ART has cd4 350+

art <- 249 # Double check !!
cb_art_y1 <- 310.84
cb_art <- 245.42


# Total costs: 

# Scenario 1

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
  mutate(ART_cost = mean2*art*(70/73) + mean2*cb_art*(3/73) + mean2*hosp_art, # ART cost plus hospitalization cost
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

