# Generate COST over time per capita for each Scenario
# MSahu
# July 28, 2020

library(dplyr)


# TO DO:
# 1. Inflate any costs
# 2. 


# LITERATURE COSTS

cost_testing_pc <- 8.58 # Sharma et al, Home HTC (need to inflate from 2012 USD to 2020)

hosp_200 <- 121
hosp_200_350 <- 58
hosp_350 <- 39
hosp_art <- 45  # This assumes 

art <- 682

# Scenario 1 cost per person



# ART cost


# Hospitalization cost


# Scenario 2 cost per person

# Generate testing costs

cost_testing <- data.frame("year"=2017:2060) %>% 
  mutate(cost=ifelse(year%%5==0,cost_testing_pc,0)) %>% 
  mutate(cost_cum=cumsum(cost))