library(tidyverse)
library(janitor)
library(ggplot2)
library(lubridate)

setwd("Christine Documents/DRIVE Model/SACEA/") # *******SET ME***********

# Translating Matlab indexes for compartments into R factors
hivCateg = data.frame("num" = seq(1,7,1),  
                         "hivCateg" = c("hiv.negative", "hiv.acute", "hiv.500", "hiv.350.500", "hiv.200.350", "hiv.200", "hiv.art"))
ageCateg = data.frame("num" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17), 
                     "ageCateg" = c("age.0.4", "age.5.9", "age.10.14", "age.15.19", "age.20.24", "age.25.29", "age.30.34", "age.35.39", 
                                      "age.40.44", "age.45.49", "age.50.54", "age.55.59", "age.60.64", "age.65.69", "age.70.74", "age.75.79", "all.ages"))
hpvCcCateg = data.frame("num" = seq(1,10,1), 
                        "hpvCcCateg" = c("hpv.susceptible", "hpv.infected", "cin1", "cin2", "cin3", "local.cc", "regional.cc", "distant.cc", "hysterectomy", "hpv.immune"))
vaxScreenCateg = data.frame("num" = seq(1,4,1), 
                            "vaxScreenCateg" = c("neither", "vax", "screen", "both"))
paramCateg = data.frame("num" = seq(1,25,1), 
                        "paramVal" = c('6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , '6_12' , '6_13' , '6_15' , '6_20' , 
                                       '6_21' , '6_22' , '6_26' , '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' ,  
                                       '6_41' , '6_42' , '6_45' , '6_47'))

deaths = read.csv("deaths.csv")
screenTreat = read.csv("screenTreat.csv")
hpvHealthState = read.csv("hpvHealthState.csv")
hivHealthState = read.csv("hivHealthState.csv")
totalPerAge = read.csv("totalPerAge.csv")
vax = read.csv("vax.csv")

# DEATHS 

deaths_clean <- deaths %>% 
  left_join(., ageCateg, by=c("age" = "num")) %>% 
  left_join(., paramCateg, by=c("paramNum" = "num")) %>% 
  left_join(., (totalPerAge %>% rename(N=count)), by=c("sceNum", "paramNum", "year", "age")) %>% 
  mutate(deathCateg = case_when(categ == 1 ~ "ccDeath", 
                                categ == 2 ~ "hivDeath", 
                                categ == 3 ~ "allDeath"), 
         sceNum = sceNum - 1) %>% # the scenario actually starts at 0, but matlab indices start at 1
  select(sceNum, paramVal, year, ageCateg, N, deathCateg, count) %>% unique %>% 
  spread(., deathCateg, count) 

# SCREENING AND TREATMENT

screenTreat_clean <- screenTreat %>% 
  left_join(., ageCateg, by=c("age" = "num")) %>% 
  left_join(., paramCateg, by=c("paramNum" = "num")) %>% 
  mutate(screenTreatCateg = case_when(categ == 1 ~ "numScreen", 
                                      categ == 2 ~ "numLEEP", 
                                      categ == 3 ~ "numCryo", 
                                      categ == 4 ~ "numThrml", 
                                      categ == 5 ~ "numHyst"), 
         sceNum = sceNum - 1) %>% 
  select(sceNum, paramVal, year, ageCateg, screenTreatCateg, count) %>% unique %>% 
  spread(., screenTreatCateg, count)

# HPV HEALTH STATES

hpvHealthState_clean <- hpvHealthState %>% 
  left_join(., ageCateg, by=c("age" = "num")) %>% 
  left_join(., paramCateg, by=c("paramNum" = "num")) %>% 
  left_join(., hpvCcCateg, by=c("categ" = "num")) %>% 
  mutate(sceNum = sceNum - 1) %>% 
  select(sceNum, paramVal, year, ageCateg, hpvCcCateg, count) %>% unique %>% 
  spread(., hpvCcCateg, count)

# HIV HEALTH STATES 

hivHealthState_clean <- hivHealthState %>% 
  left_join(., ageCateg, by=c("age" = "num")) %>% 
  left_join(., paramCateg, by=c("paramNum" = "num")) %>% 
  left_join(., hivCateg, by=c("categ" = "num")) %>% 
  mutate(sceNum = sceNum - 1) %>% 
  select(sceNum, paramVal, year, ageCateg, hivCateg, count) %>% unique %>% 
  spread(., hivCateg, count)

# VACCINATIONS
vax_clean <- vax %>% 
  left_join(., ageCateg, by=c("age" = "num")) %>% 
  left_join(., paramCateg, by=c("paramNum" = "num")) %>% 
  mutate(sceNum = sceNum - 1, 
         bivalCount = case_when(year < 2021 ~ count, TRUE ~ 0), # before 2021, all counts are bivalent
         nonavalCount = case_when(year >= 2021 ~ count, TRUE ~ 0)) %>% # after 2021, all counts are nonavalent
  select(sceNum, paramVal, year, ageCateg, bivalCount, nonavalCount) %>% unique 

# Combine all the dfs
combined <- deaths_clean %>% 
  left_join(., screenTreat_clean, by=c("sceNum", "paramVal", "year", "ageCateg")) %>% 
  left_join(., hpvHealthState_clean, by=c("sceNum", "paramVal", "year", "ageCateg")) %>% 
  left_join(., hivHealthState_clean, by=c("sceNum", "paramVal", "year", "ageCateg")) %>% 
  left_join(., vax_clean, by=c("sceNum", "paramVal", "year", "ageCateg")) %>% 
  mutate(numScreen = case_when(is.na(numScreen) ~ 0, TRUE ~ numScreen), # when these values are NA, convert to 0
         numLEEP = case_when(is.na(numLEEP) ~ 0, TRUE ~ numLEEP), 
         numCryo = case_when(is.na(numCryo) ~ 0, TRUE ~ numCryo), 
         numThrml = case_when(is.na(numThrml) ~ 0, TRUE ~ numThrml), 
         numHyst = case_when(is.na(numHyst) ~ 0, TRUE ~ numHyst)) %>% 
  select(sceNum, paramVal, year, ageCateg, N, hpv.susceptible, hpv.immune, hpv.infected, cin1, cin2, cin3, 
          local.cc, regional.cc, distant.cc, hysterectomy, hiv.negative, hiv.art, hiv.500, hiv.350.500,
          hiv.200.350, hiv.200, hiv.acute, nonavalCount, bivalCount, numScreen, numLEEP, numCryo, numThrml, numHyst,
          allDeath, ccDeath, hivDeath) %>% unique 

# Loop through scenarios and output the final files per scenario 

for (i in seq(0,9,1)) {
  sceDf <- combined %>% 
    filter(sceNum == i)
  
  ifelse(!dir.exists(paste0(getwd(), "/Outputs")), dir.create(paste0(getwd(), "/Outputs")), FALSE)
  
  write.csv(sceDf, paste0(getwd(), "/Outputs/modelResultsForCea_S", i, ".csv"), row.names = FALSE)
}
