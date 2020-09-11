
# Compile cases averted across DO-ART scenarios, combined for Males and Females
# MSahu
# September 10, 2020

# NOTE that annual incidence is calculated per 100 HIV-negative persons 


##############################################################################

# Setup

source("00_master.R")
library("dplyr")

##############################################################################

# IMPORT POPULATION (males and females combined, age 15-79)

##############################################################################

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
rm(pop)

##############################################################################

# IMPORT INCIDENCE (males and females combined, age 15-79)

##############################################################################

for (x in 1:3) {

incidence <- read.csv(paste0(main_path,"Scenario",x,"/HIV_incidence_combined_aged15-79.csv"), header=F) %>% 
  setNames(paste0("s",-3:25)) %>% 
  rename(year=1,
         mean=2,
         min=3,
         max=4)   %>% 
  filter(year>=2020) 

assign(paste0("incidence",x),incidence)

}

##############################################################################

# IMPORT PREVALENCE (males and females combined, age 15-79)

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


##############################################################################

# REVISED METHOD: Calculate raw and percentage cases averted ; discounted

##############################################################################

# HIV-susceptible population is 1 - Prevalence


for (x in 1:3) {
  
  hiv_neg <- ((1- (get(paste0("prev_ART_scen",x))[,-1] +    # HIV-negative is 1 - HIV Prevalence
                     get(paste0("prev_200_scen",x))[,-1] +
                     get(paste0("prev_200_350_scen",x))[,-1] +
                     get(paste0("prev_350_scen",x))[,-1] )) 
                *(get(paste0("pop_scen",x))[,-1]) ) %>%   # multiply by population
    mutate(year=2020:2060) %>% 
    select(year, everything()) 
  
  assign(paste0("hiv_neg_scen",x),hiv_neg)
  
}

# New cases

for (x in 1:3) {
  
  cases <- ((get(paste0("incidence",x))[,-1]*get(paste0("hiv_neg_scen",x))[,-1])/100)  %>% # Get cases given incidence rate per 100   
    
    # DISCOUNT
    mutate(year_discount=0:40,   # set 2020 to Year 0
               discount_amt=discount(discount_rate,year_discount)) %>% 
    transmute_at(vars(c(mean,min,max,4:28)),~discounter(.,discount_amt)) %>% 
     
    mutate(year=2020:2060) %>% 
    select(year, everything()) 
     
  assign(paste0("cases_scen",x),cases)
  
}

# Cumulative cases

for (x in 1:3) {
  
  cases_cum <- cumsum(get(paste0("cases_scen",x))[,-1])  %>% 
    mutate(year=2020:2060) %>% 
    select(year, everything()) 
  
  assign(paste0("cases_cum_scen",x),cases_cum)
  
}

###################################################################################################################

# Annual raw cases averted 

cases_averted_raw_scen2 <- (cases_scen1[,-1] - cases_scen2[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

cases_averted_raw_scen3 <- (cases_scen1[,-1] - cases_scen3[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Annual percent cases averted

cases_averted_pct_scen2 <- (((cases_scen1[,-1] - cases_scen2[,-1] )/cases_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

cases_averted_pct_scen3 <- (((cases_scen1[,-1] - cases_scen3[,-1] )/cases_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative raw cases averted 

cases_cum_averted_raw_scen2 <- ( cases_cum_scen1[,-1] - cases_cum_scen2[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

cases_cum_averted_raw_scen3 <- ( cases_cum_scen1[,-1] - cases_cum_scen3[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative percent cases averted

cases_cum_averted_pct_scen2 <- (((cases_cum_scen1[,-1] - cases_cum_scen2[,-1] )/cases_cum_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

cases_cum_averted_pct_scen3 <- (((cases_cum_scen1[,-1] - cases_cum_scen3[,-1] )/cases_cum_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

###################################################################################################################

# Export CSVs

csv_list <- list("cases_scen1",
                 "cases_scen2",
                 "cases_scen3",
                 "cases_cum_scen1",
                 "cases_cum_scen2",
                 "cases_cum_scen3",
                 "cases_averted_raw_scen2",
                 "cases_averted_raw_scen3",
                 "cases_averted_pct_scen2",
                 "cases_averted_pct_scen3",
                 "cases_cum_averted_raw_scen2",
                 "cases_cum_averted_raw_scen3",
                 "cases_cum_averted_pct_scen2",
                 "cases_cum_averted_pct_scen3")


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"effects/cases/",x,".csv")))

# Remove excess DFs

rm(list=ls())
