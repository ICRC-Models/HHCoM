# Calculate QALYs
# MSahu
# September 9, 2020

# QALY Calculation: Years of Life * Utility Value
# This script calculates average annual QALYs gained between interventions

# TO DO:
# Check if CD4+ and ART status of incident cases for model is available

library(dplyr)

# Set up utility weight parameters

wt_hiv_neg <- 1
wt_cd4_350 <- 0.94
wt_cd4_250_350 <- 0.82
wt_cd4_200 <- 0.7
wt_art <- 0.94
wt_dead <- 0

# Current strategy:
#1. Number of prevalent cases * utility weight (by CD4+ count), per capita
#2. Number of incident cases * utility weight (Assuming CD4+ >350 for all), per capita
#3. Number of HIV negative cases * utility weight  [can be calculated as population minus inc minus prev minus deaths], per capita
#4. Sum utilities
#5. Discount

# Then Ask Paul Revill - because actually should be cumulative ! (or are births available?)
# Also check if reasonable to assume all incident cases have CD4+ >350 


##############################################################################

# IMPORT PREVALENCE (males and females combined), by CD4+ count and ART status

##############################################################################

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



# Mortality by year, for hypothetical population  (Mortality rates already given per 100k)

for (x in 1:3) {
  
  mortality <- read.csv(paste0(main_path,"Scenario",x,"/HIV_mortality_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) %>% 
    filter(year>=2020) %>% 
    select(-year)
  
  assign(paste0("mortality",x),mortality)
  
}

rm(mortality)


# QALYs by year: Cumulative and discounted



