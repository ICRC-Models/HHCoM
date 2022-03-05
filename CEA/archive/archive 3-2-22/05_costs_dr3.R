######################################################################

# Calculate average annual costs per year
# Author: Mita Sahu
# Date: September 8, 2020

# NOTE: Must be connected to IHME VPN for currency conversion function 

# TO do: 
# 1. Diagnosed / Undiagnosed for Scenario 3? (asked Cara)
# 2. Add Scenario 2a

######################################################################

# Setup ; Read inflation converter from H drive

setwd("/")

library(pacman)

source('H:/repos/fgh/FGH_2019/04_functions/currency_conversion.R')

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/")

#######################################################################

# Inflate literature costs ; assign parameter values

#######################################################################

cost_params <- read.csv(paste0(helper_path,"cost_params.csv")) %>% 
  mutate(iso3="USA") 

cost_inflated <- currency_conversion(data = cost_params,
                      col.value = "param_val",
                      col.loc = "iso3", 
                      base.year = 2020, 
                      col.currency="currency", 
                      col.currency.year="year")

costs <- cost_inflated[["param_val"]]
names(costs) <- cost_inflated[["param_name"]]

# MAXIMUM ALLOWABLE COSTS 
# costs[["cb_art_y1"]] = 2800
# costs[["cb_art_sub"]] = 2800

########################################################################

# IMPORT POPULATION

########################################################################

# Full population

for (x in 1:3) {
  
  pop <- read.csv(paste0(main_path,"Scenario",x,"/PopulationSize_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("pop_scen",x),pop)
  
}

# Male population  

for (x in 1:3) {
  
  pop <- read.csv(paste0(main_path,"Scenario",x,"/PopulationSize_males_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Mpop_scen",x),pop)
  
}

# Female population  

for (x in 1:3) {
  
  pop <- read.csv(paste0(main_path,"Scenario",x,"/PopulationSize_females_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("Fpop_scen",x),pop)
  
}

########################################################################

# IMPORT PREVALENCE

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

# Prevalence on ART, males

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

# Prevalence on ART, females

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

########################################################################

# IMPORT INCIDENCE and calculate NEW CASES (PER CAPITA)

# Incidence rates are given per 100 HIV-negative persons
# We need new cases per total population, gender-specific

# Multiply incidence per HIV-susceptible * proportion HIV-negative

########################################################################

# Incidence rates per capita among HIV-susceptibles

sex <- c("males","females")

for (i in sex) {

  for (x in 1:3) {
  
  incidence <- as.data.frame(read.csv(paste0(main_path,"Scenario",x,"/HIV_incidence_",i,"_aged15-79.csv"), header=F)[,-1] / 100 ) %>% 
    setNames(paste0("s",-2:25)) %>% 
    mutate(year=1925:2060) %>%
    filter(year>=2020) %>%
    select(year, everything()) %>%
    rename(mean=2,
           min=3,
           max=4)  
  
  assign(paste0("incidence_",i,x),incidence)
  }
}

# Calculate proportion HIV-negative as 1- Prevalence, by sex and aggregate

for (i in sex) {
  
  for (x in 1:3) {
    
    hiv_neg <- as.data.frame(1-read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_",i,"_aged15-79_all_HIV.csv"), header=F)[,-1]) %>% 
      setNames(paste0("s",-2:25)) %>% 
      mutate(year=1925:2060) %>%
      filter(year>=2020) %>%
      select(year, everything()) %>%
      rename(mean=2,
             min=3,
             max=4)  
    
    assign(paste0("hiv_neg_",i,x),hiv_neg)
  }
}

for (x in 1:3) {
  
  hiv_neg <- as.data.frame(1-read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_all_HIV.csv"), header=F)[,-1]) %>% 
    setNames(paste0("s",-2:25)) %>% 
    mutate(year=1925:2060) %>%
    filter(year>=2020) %>%
    select(year, everything()) %>%
    rename(mean=2,
           min=3,
           max=4)  
  
  assign(paste0("hiv_neg_",x),hiv_neg)
}

# Calculate new cases per total population of M/F by multiplying incidence rates among HIV-negative by proportion HIV-negative


for (i in sex) {
  
  for (x in 1:3) {
    
    new_cases <- as.data.frame(get(paste0("hiv_neg_",i,x))[,-1] * get(paste0("incidence_",i,x))[,-1]) %>% 
      addYearCol(., horizon_year = horizon_year) 
    
    assign(paste0("cases_",i,x),new_cases)
  }
}


############################################################################################

# Remove excess DFs

rm(prev_ART, prev_200, prev_200_350, prev_350, prev_ART_males, prev_ART_females, incidence, hiv_neg)


############################################################################################

# HOSPITALIZATION COSTS

# Multiply HIV prevalence by average annual cost of hospitalization, by CD4+ count

############################################################################################

# All 3 Scenarios follow same formula:

for (x in 1:3) {

hosp <- as.data.frame(
              (get(paste0("prev_200_scen",x))[,-1]) * costs["hosp_200"] +
              (get(paste0("prev_200_350_scen",x))[,-1]) * costs["hosp_200_350"] +
              (get(paste0("prev_350_scen",x))[,-1]) * costs["hosp_350"] +
              (get(paste0("prev_ART_scen",x))[,-1]) * costs["hosp_art"] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

assign(paste0("hosp_scen",x),hosp)

}

rm(hosp)


####################################################################################################

# TESTING COSTS (PER CAPITA)

# For Scenario 2 and 3 only, assume one testing campaign every 5 years
# This testing campaign covers HIV negatives plus undiagnosed HIV positives
# Assumes that 90% of people will be reached.

# INPUT: Cara's files on numbers of people tested "Raw_HTC" ; cost and coverage assumptions
# OUTPUT: Total testing cost per year, undiscounted, by scenario

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

# HIV-negatives
test_scen2_hivNeg <- read.csv(paste0(main_path,"Scenario2/Raw_HTC_combined_hivNeg.csv"), 
                              header = F)[,-1] * costs["home_testing_pc_neg"] # MULTIPLY BY cost of negative test 

# HIV-positive undiagnosed
test_scen2_hivPos <- read.csv(paste0(main_path,"Scenario2/Raw_HTC_combined_hivUndiag.csv"), 
                              header = F)[,-1] * costs["home_testing_pc_pos"] # MULTIPLY BY cost of positive test 

test_scen2 <- ( (test_scen2_hivNeg + test_scen2_hivPos) / (pop_scen2[,-1]) ) %>%  # PER CAPITA 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1, min=2, max=3) %>% 
  mutate(year=row_number()+2019) %>% 
  select(year,everything())

rm(test_scen2_hivNeg, test_scen2_hivPos)
  

# Scenario 3 includes a testing campaign every 5 years, per person, starting 2020

# ?? Follow up with Cara - is this done for Scenario 3
# PLACEHOLDER: all zeros

test_scen3 <- as.data.frame(matrix(0, ncol = 28, nrow = 41)) %>% 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1,
                min=2,
                max=3) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything())

####################################################################################################

# ART COSTS (PER CAPITA)

# Costs for prevalent cases + costs for incident cases (for 1/2 year)

# Multiply prevalence of people on ART by the ART costs, using scalar for percentage on community 
# versus clinic-based ART (disaggregated for males and females)

# Incident cases are assumed to be enrolled using the model assumption for % cases starting ART
# Costs are for 1/2 year (disaggregated by M/F)

####################################################################################################

# Generate cb-ART costs for the first year and subsequent years 
                                   
cbART_costs <- data.frame(matrix(costs["cb_art_y1"],ncol = 28, nrow = 1)) %>% 
  rbind(data.frame(matrix(costs["cb_art_sub"],ncol = 28, nrow = 40))) %>% 
  setNames(paste0("s",-2:25)) %>%
  dplyr::rename(mean=1,
                min=2,
                max=3) %>%
  mutate(year=2020:2060) %>% 
  select(year,everything())

######################################
# SCENARIO 1 : Standard of care ONLY #
######################################

# MALES

mART_scen1 <- as.data.frame(   
  # Prevalent cases * standard of care costs   
  (Mprev_ART_scen1[,-1] * costs["standard_art"]) + 
  # Incident cases * % initiating ART * 1/2 year * standard of care costs
  (cases_males1[,-1] * males_ART_scen1 * costs["standard_art"] * 0.5 ))  

# FEMALES

fART_scen1 <- as.data.frame(   
  # Prevalent cases * standard of care costs 
  (Fprev_ART_scen1[,-1] * costs["standard_art"]) + 
  # Incident cases * % initiating ART * 1/2 year * standard of care costs
  (cases_females1[,-1] * females_ART_scen1 * costs["standard_art"] * 0.5 ))  

##########################################
# SCENARIO 2 : Standard of care + cb-ART #
##########################################

# MALES 

mART_scen2 <- as.data.frame(
  # Prevalent cases * enrollment scalar * ART costs 
  (Mprev_ART_scen2[,-1] * males_enrolment_clinic * costs["standard_art"]) + # Standard of care cost
  (Mprev_ART_scen2[,-1] * males_enrolment_cbART * cbART_costs[,-1]) +   # Community-based ART cost
  # Incident cases * % initiating ART * enrollment scalar * 1/2 year * ART costs
  (cases_males2[,-1] * males_ART_scen2 * males_enrolment_clinic * costs["standard_art"] * 0.5 ) + # Standard of care cost
  (cases_males2[,-1] * males_ART_scen2 * males_enrolment_cbART * cbART_costs[,-1] * 0.5 )) # Community-based ART cost

# FEMALES

fART_scen2 <- as.data.frame(
  # Prevalent cases * enrollment scalar * ART costs 
  (Fprev_ART_scen2[,-1] * females_enrolment_clinic * costs["standard_art"]) + # Standard of care cost
  (Fprev_ART_scen2[,-1] * females_enrolment_cbART * cbART_costs[,-1]) +   # Community-based ART cost
  # Incident cases * % initiating ART * enrollment scalar * 1/2 year * ART costs
  (cases_females2[,-1] * females_ART_scen2 * females_enrolment_clinic * costs["standard_art"] * 0.5 ) + # Standard of care cost
  (cases_females2[,-1] * females_ART_scen2 * females_enrolment_cbART * cbART_costs[,-1] * 0.5 )) # Community-based ART cost

######################################################
# SCENARIO 3 (base case) : Standard of care + cb-ART #
######################################################

# MALES 

mART_scen3 <- as.data.frame(
  # Prevalent cases * enrollment scalar * ART costs 
  (Mprev_ART_scen3[,-1] * males_enrolment_clinic * costs["standard_art"]) + # Standard of care cost
  (Mprev_ART_scen3[,-1] * males_enrolment_cbART * cbART_costs[,-1]) +   # Community-based ART cost
  # Incident cases * % initiating ART * enrollment scalar * 1/2 year * ART costs
  (cases_males3[,-1] * males_ART_scen3 * males_enrolment_clinic * costs["standard_art"] * 0.5 ) + # Standard of care cost
  (cases_males3[,-1] * males_ART_scen3 * males_enrolment_cbART * cbART_costs[,-1] * 0.5 )) # Community-based ART cost

# FEMALES

fART_scen3 <- as.data.frame(
  # Prevalent cases * enrollment scalar * ART costs 
  (Fprev_ART_scen3[,-1] * females_enrolment_clinic * costs["standard_art"]) + # Standard of care cost
  (Fprev_ART_scen3[,-1] * females_enrolment_cbART * cbART_costs[,-1]) +   # Community-based ART cost
  # Incident cases * % initiating ART * enrollment scalar * 1/2 year * ART costs
  (cases_females3[,-1] * females_ART_scen3 * females_enrolment_clinic * costs["standard_art"] * 0.5 ) + # Standard of care cost
  (cases_females3[,-1] * females_ART_scen3 * females_enrolment_cbART * cbART_costs[,-1] * 0.5 )) # Community-based ART cost

############################  
# COMBINED MALE AND FEMALE #
############################

# Generate weight - proportion male

for (x in 1:3) {
  
  M_wt <- ( (get(paste0("Mpop_scen",x))[,-1]) / (get(paste0("pop_scen",x))[,-1]))%>% 
    mutate(year=2020:2060) %>% 
    select(year,everything())
  
  assign(paste0("Mwt_scen",x),pop)
  
}

# Combine male and female, by applying weight for males and females

for (x in 1:3) {

ART <- (( get(paste0("mART_scen",x)) * (M_wt[,-1]) )  +       # Male cost * male wt
       (  get(paste0("fART_scen",x)) * ((1-M_wt)[,-1]) )) %>% # Female cost * (1- male wt)
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

assign(paste0("ART_scen",x),ART)

}

# Remove excess

rm(pop, ART, M_wt, mART_scen1,fART_scen1, mART_scen2, fART_scen2, mART_scen3, fART_scen3, Mwt_scen1, Mwt_scen2, Mwt_scen3)


####################################################################################################

# TOTAL COSTS

# Sum the Hospitalization, Testing, and ART Costs
# Then discount

####################################################################################################

# 1. Annual per capita (discounted)
# 2. Annual total (discounted)
# 3. Cumulative per capita (discounted)
# 4. Cumulative total (discounted)

# 5. Incremental costs (annual per capita)
# 6. Incremental costs  (annual total)
# 7. Incremental costs  (cumulative per capita)
# 8. Incremental costs  (cumulative total)


# 1. ANNUAL PER CAPITA (discounted)

for (x in 1:3) {

cost <- (get(paste0("hosp_scen",x))[,-1] + 
           get(paste0("test_scen",x))[,-1]+
           get(paste0("ART_scen",x))[,-1]) %>% 
  
  # DISCOUNT
  discount(., discount_rate = discount_rate) %>%  
  
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max
  
assign(paste0("cost_annual_pc_scen",x),cost)

}

# 2. ANNUAL TOTAL (discounted)

for (x in 1:3) {
  
  cost <- ( get(paste0("cost_annual_pc_scen",x))[,-1] * get(paste0("pop_scen",x))[,-1] ) %>% 
    addYearCol(., horizon_year = horizon_year) %>% 
    recalcFuns(.) # recalculate mean, min, max
  
  assign(paste0("cost_annual_total_scen",x),cost)
  
}


# 3. CUMULATIVE, per capita (discounted)

for (x in 1:3) {
  
  cost <- (  get(paste0("hosp_scen", x))[,-1] + 
             get(paste0("test_scen", x))[,-1] +
             get(paste0("ART_scen", x))[,-1] ) %>% 
    
    # DISCOUNT
    discount(., discount_rate = discount_rate) %>%  
    
    # Cumulative
    transmute_all(~cumsum(.)) %>% 
    
    addYearCol(., horizon_year = horizon_year) %>% 
    recalcFuns(.) # recalculate mean, min, max
    
  assign(paste0("cost_cum_pc_scen",x),cost)
  
}


# 4. CUMULATIVE, total (discounted)

for (x in 1:3) {
  
  cost <- (get(paste0("cost_cum_pc_scen",x))[,-1]*get(paste0("pop_scen",x))[,-1]) %>% 
    addYearCol(., horizon_year = horizon_year) %>% 
    recalcFuns(.) # recalculate mean, min, max
  
  assign(paste0("cost_cum_total_scen",x),cost)
  
}

# 5. Incremental costs (annual per capita)

inc_costs_annual_pc_scen2 <- ( cost_annual_pc_scen2[,-1] - cost_annual_pc_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max 

inc_costs_annual_pc_scen3 <- ( cost_annual_pc_scen3[,-1] - cost_annual_pc_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# 6. Incremental costs  (annual total)

inc_costs_annual_total_scen2 <- ( cost_annual_total_scen2[,-1] - cost_annual_total_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

inc_costs_annual_total_scen3 <- ( cost_annual_total_scen3[,-1] - cost_annual_total_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# 7. Incremental costs  (cumulative per capita)

inc_costs_cum_pc_scen2 <- ( cost_cum_pc_scen2[,-1] - cost_cum_pc_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

inc_costs_cum_pc_scen3 <- ( cost_cum_pc_scen3[,-1] - cost_cum_pc_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# 8. Incremental costs  (cumulative total)

inc_costs_cum_total_scen2 <- ( cost_cum_total_scen2[,-1] - cost_cum_total_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

inc_costs_cum_total_scen3 <- ( cost_cum_total_scen3[,-1] - cost_cum_total_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

##############################################################################################

# Export CSVs

csv_list <- list("hosp_scen1",
                 "hosp_scen2",
                 "hosp_scen3",
                 "ART_scen1",
                 "ART_scen2",
                 "ART_scen3",
                 "test_scen1",
                 "test_scen2",
                 "test_scen3",
                 "cost_annual_pc_scen1",
                 "cost_annual_pc_scen2",
                 "cost_annual_pc_scen3",
                 "cost_annual_total_scen1",
                 "cost_annual_total_scen2",
                 "cost_annual_total_scen3",
                 "cost_cum_pc_scen1",
                 "cost_cum_pc_scen2",
                 "cost_cum_pc_scen3",
                 "cost_cum_total_scen1",
                 "cost_cum_total_scen2",
                 "cost_cum_total_scen3",
                 "inc_costs_annual_pc_scen2",
                 "inc_costs_annual_pc_scen3",
                 "inc_costs_annual_total_scen2",
                 "inc_costs_annual_total_scen3",
                 "inc_costs_cum_pc_scen2",
                 "inc_costs_cum_pc_scen3",
                 "inc_costs_cum_total_scen2",
                 "inc_costs_cum_total_scen3"
           )


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"results/dr3/costs/",x,".csv"), row.names = F ))

# Remove excess DFs, source for next script

rm(list=ls()[! ls() %in% c("main_path","cea_path","helper_path")])
source(paste0(helper_path, "helper.R"))

