######################################################################

# Calculate average annual costs per year
# Author: Mita Sahu
# Date: September 8, 2020

# NOTE: Must be connected to IHME VPN for currency conversion function 


# !!!! NEED TO RESOLVE SCENARIO 3
# !!!! Need to add ART and hospitalization costs for incident cases

######################################################################

source("00_master.R")
library(pacman)


source('H:/repos/fgh/FGH_2019/04_functions/currency_conversion.R')

library(dplyr)

# TO DO:
# 1. Double check if reasonable to assume everyone on ART is CD4>350
# 2. Resolve scenario 3
# 3. Calculate QALYS!!
# 4. Add costs of incident cases.

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
    dplyr::rename(year=1,
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
    dplyr::rename(year=1,
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
    dplyr::rename(year=1,
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
    dplyr::rename(year=1,
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
              prev_350_scen1[,-1] * costs["hosp_350"] +
              prev_ART_scen1[,-1] * costs["hosp_art"] ) %>% 
      mutate(year=2020:2060) %>% 
      select(year,everything())

# Scenario 2

hosp_scen2 <- as.data.frame(
  prev_200_scen2[,-1] * costs["hosp_200"] +
    prev_200_350_scen2[,-1] * costs["hosp_200_350"] +
    prev_350_scen2[,-1] * costs["hosp_350"] +
    prev_ART_scen2[,-1] * costs["hosp_art"] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 3

hosp_scen3 <- as.data.frame(
  prev_200_scen3[,-1] * costs["hosp_200"] +
    prev_200_350_scen3[,-1] * costs["hosp_200_350"] +
    prev_350_scen3[,-1] * costs["hosp_350"] +
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
  dplyr::rename(mean=1,
         min=2,
         max=3) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 2 includes a testing campaign every 5 years, per person, starting 2020

test_scen2 <- as.data.frame(matrix(rep(c(costs["home_testing_pc"],0,0,0,0)), ncol = 28, nrow = 40)) %>% 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1,
         min=2,
         max=3) %>% 
  rbind(costs["home_testing_pc"]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Scenario 3 includes a testing campaign every 5 years, per person, starting 2020

test_scen3 <- as.data.frame(matrix(rep(c(costs["home_testing_pc"],0,0,0,0)), ncol = 28, nrow = 40)) %>% 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1,
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
    dplyr::rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Mprev_ART_scen",x),prev_ART_males)
  
}

for (x in 1:3) {
  
  prev_ART_females<- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_females_aged15-79_onART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Fprev_ART_scen",x),prev_ART_females)
  
}

rm(prev_ART_males,prev_ART_females)


# Generate AVERAGE cb-ART costs for the first year and costs for subsequent years 

cbART_costs_pc <- data.frame(year=2020:2060) %>% 
  mutate(cb_ART_cost= 
           (costs["cb_art_y1"] + (year-2020)*costs["cb_art_sub"] )/  # Numerator is 1st year costs, plus all subsequent years  
           (year-2020+1))   # Denominator is total number of years
                                   
cbART_costs <- data.frame(matrix(cbART_costs_pc[,2],ncol = 28, nrow = 41)) %>% 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1,
                min=2,
                max=3) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())
  
rm(cbART_costs_pc)


# Multiply based on percentages enrolled in clinic versus community-based ART (defined in master script)

# Males

mART_scen1 <- as.data.frame(
  (Mprev_ART_scen1[,-1] * males_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Mprev_ART_scen1[,-1] * males_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

mART_scen2 <- as.data.frame(
  (Mprev_ART_scen2[,-1] * males_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Mprev_ART_scen2[,-1] * males_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

# NEED TO RESOLVE SCENARIO 3

mART_scen3 <- as.data.frame(
  (Mprev_ART_scen3[,-1] * males_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Mprev_ART_scen3[,-1] * males_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

# Females

fART_scen1 <- as.data.frame(
  (Fprev_ART_scen1[,-1] * females_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Fprev_ART_scen1[,-1] * females_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

fART_scen2 <- as.data.frame(
  (Fprev_ART_scen2[,-1] * females_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Fprev_ART_scen2[,-1] * females_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

# NEED TO RESOLVE SCENARIO 3

fART_scen3 <- as.data.frame(
  (Fprev_ART_scen3[,-1] * females_enrolment_clinic * costs["standard_art"]) + #Standard of care cost
  (Fprev_ART_scen3[,-1] * females_enrolment_cbART * cbART_costs[,-1]))   # Community-based ART cost

# Combined

ART_scen1 <- (mART_scen1 + fART_scen1) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

ART_scen2 <- (mART_scen2 + fART_scen2) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# NEED TO RESOLVE SCENARIO 3

ART_scen3 <- (mART_scen3 + fART_scen3) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

# Remove excess

rm(mART_scen1,fART_scen1, mART_scen2, fART_scen2, mART_scen3, fART_scen3)


####################################################################################################

# TOTAL COSTS

# Sum the Hospitalization, Testing, and ART Costs
# Then discount

####################################################################################################


# Scenario 1

COST_scen1 <- (hosp_scen1[,-1] + test_scen1[,-1] + ART_scen1[,-1]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) %>% 
  # DISCOUNT
  mutate(year_discount=0:40,   # set 2020 to Year 0
         discount_amt=discount(discount_rate,year_discount)) %>% 
  mutate_at(vars(c(mean,min,max,5:29)),~discounter(.,discount_amt)) %>% 
  select(-year_discount,-discount_amt)

write.csv(COST_scen1, paste0(cea_path,"costs/costs_scen1.csv"))

# Scenario 2

COST_scen2 <- (hosp_scen2[,-1] + test_scen2[,-1] + ART_scen2[,-1]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) %>% 
  # DISCOUNT
  mutate(year_discount=0:40,   # set 2020 to Year 0
         discount_amt=discount(discount_rate,year_discount)) %>% 
  mutate_at(vars(c(mean,min,max,5:29)),~discounter(.,discount_amt)) %>% 
  select(-year_discount,-discount_amt)

write.csv(COST_scen2, paste0(cea_path,"costs/costs_scen2.csv"))

# Scenario 3

COST_scen3 <- (hosp_scen3[,-1] + test_scen3[,-1] + ART_scen3[,-1]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) %>% 
  # DISCOUNT
  mutate(year_discount=0:40,   # set 2020 to Year 0
         discount_amt=discount(discount_rate,year_discount)) %>% 
  mutate_at(vars(c(mean,min,max,5:29)),~discounter(.,discount_amt)) %>% 
  select(-year_discount,-discount_amt)

write.csv(COST_scen3, paste0(cea_path,"costs/costs_scen3.csv"))
