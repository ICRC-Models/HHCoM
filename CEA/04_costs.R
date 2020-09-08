######################################################################

# Calculate average annual costs at 2030, 2040, 2050, and 2060
# Author: Mita Sahu
# Date: September 8, 2020

# NOTE: Must be connected to IHME VPN for currency conversion function 

######################################################################

library(dplyr)
library(pacman)
source("00_master.R")
source('H:/repos/fgh/FGH_2019/04_functions/currency_conversion.R')


# TO DO:
# 1. Inflate any costs -- DONE
# 2. Discount
# 3. Double check if reasonable to assume everyone on ART is CD4>350
# 4. Check ART costs
# 5. Include startup cost cb_art ,
# 6. Calculate QALYS!!

#######################################################################

# Inflate literature costs ; assign parameter values

#######################################################################

cost_params <- read.csv(paste0(cea_path,"costs/cost_params.csv")) %>% 
  mutate(iso3="USA") 

cost_inflated <- currency_conversion(data = cost_params,
                      col.value = "param_val",
                      col.loc = "iso3", 
                      base.year = 2020, 
                      col.currency="currency", 
                      col.currency.year="year")

costs <- cost_inflated[["param_val"]]
names(costs) <- cost_inflated[["param_name"]]

########################################################################

# IMPORT PREVALENCE (males and females combined)

########################################################################

# Prevalence on ART

for (x in 1:3) {
  
  prev_ART<- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_onART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_ART_scen",x),prev_ART)
  
}

# Prevalence CD4 below 200, not on ART

for (x in 1:3) {
  
  prev_200 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_below200_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_200_scen",x),prev_200)
  
}

# Prevalence CD4 200-350, not on ART

for (x in 1:3) {
  
  prev_200_350 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_200-350_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_200_350_scen",x),prev_200_350)
  
}

# Prevalence CD4 350+, not on ART 

for (x in 1:3) {
  
  prev_350 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_350plus_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_350_scen",x),prev_350)
  
}

# Remove excess DFs

rm(prev_ART, prev_200, prev_200_350, prev_350)


##################################################################################

# HOSPITALIZATION COSTS

# Multiply HIV prevalence by average annual cost of hospitalization, by CD4+ count

##################################################################################

# Scenario 1

hosp_scen1 <- as.data.frame(
              prev_200_scen1[,-1] * costs["hosp_200"] +
              prev_200_350_scen1[,-1] * costs["hosp_200_350"] +
              prev_350_scen1[,-1] * costs["hosp_200_350"] +
              prev_ART_scen1[,-1] * costs["hosp_art"] ) %>% 
      mutate(year=2020:2060) %>% 
      select(year,everything())

# Scenario 2

hosp_scen2 <- as.data.frame(
  prev_200_scen2[,-1] * costs["hosp_200"] +
    prev_200_350_scen2[,-1] * costs["hosp_200_350"] +
    prev_350_scen2[,-1] * costs["hosp_200_350"] +
    prev_ART_scen2[,-1] * costs["hosp_art"] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 3

hosp_scen3 <- as.data.frame(
  prev_200_scen3[,-1] * costs["hosp_200"] +
    prev_200_350_scen3[,-1] * costs["hosp_200_350"] +
    prev_350_scen3[,-1] * costs["hosp_200_350"] +
    prev_ART_scen3[,-1] * costs["hosp_art"] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())


####################################################################################################

# TESTING COSTS

# For Scenario 2 and 3 only, assume one testing campaign every 5 years
# Note that this is an overestimate because we assume everyone (HIV positive and negative is tested)

####################################################################################################

# Scenario 1 is 0 for all

test_scen1 <- as.data.frame(matrix(0, ncol = 28, nrow = 41)) %>% 
  setNames(paste0("s",-2:25)) %>%
  rename(mean=1,
         min=2,
         max=3) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 2 includes a testing campaign every 5 years, per person, starting 2020

test_scen2 <- as.data.frame(matrix(rep(c(costs["home_testing_pc"],0,0,0,0)), ncol = 28, nrow = 40)) %>% 
  setNames(paste0("s",-2:25)) %>%
  rename(mean=1,
         min=2,
         max=3) %>% 
  rbind(costs["home_testing_pc"]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 3 includes a testing campaign every 5 years, per person, starting 2020

test_scen3 <- as.data.frame(matrix(rep(c(costs["home_testing_pc"],0,0,0,0)), ncol = 28, nrow = 40)) %>% 
  setNames(paste0("s",-2:25)) %>%
  rename(mean=1,
         min=2,
         max=3) %>% 
  rbind(costs["home_testing_pc"]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())


####################################################################################################

# ART COSTS

# Multiply prevalence of people on ART by the ART costs, using scalar for percentage on community 
# versus clinic-based ART (disaggregated for males and females)

####################################################################################################

# Import prevalence of people on ART, for males and females

for (x in 1:3) {
  
  prev_ART_males<- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_males_aged15-79_onART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Mprev_ART_scen",x),prev_ART_males)
  
}

for (x in 1:3) {
  
  prev_ART_females<- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_females_aged15-79_onART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Fprev_ART_scen",x),prev_ART_females)
  
}

rm(prev_ART_males,prev_ART_females)

# Multiply based on percentages enrolled in clinic versus community-based ART (defined in master script)

mART_scen1 <- as.data.frame(
    (Mprev_ART_scen1[,-1] * males_enrolment_clinic * costs["standard_art"]) +   # clinic ART costs
    (Mprev_ART_scen1[,-1] * males_enrolment_cbART * costs["cb_art_sub"]))
      
    prev_200_350_scen1[,-1] * costs["hosp_200_350"] +
    prev_350_scen1[,-1] * costs["hosp_200_350"] +
    prev_ART_scen1[,-1] * costs["hosp_art"] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())


#################################
# Annual costs, 2030 - Scenario 1
#################################

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

