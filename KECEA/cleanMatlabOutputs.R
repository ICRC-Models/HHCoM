# title: cleanMatlabOutputs.R
# date: 6/2/2023
# author: Christine L Hathaway
# description: FOR KENYA MODEL. Takes the CSV files spit out by Matlab (vaxCEA_multSims_mainFunction.m) and puts it into a format 
# that can be used for the economic part of the analysis. 


####### TODO
# - Make correction for the first year of vax 
# - Make a correction for 2 vs 1 dose; will depend on the scenario
# - Make a correction for the 0.7/0.9 coverage adjustment. You have to multiply by 0.9/0.7 to get the actual coverage value. 
# - Reorder the age categories so they are in increasing order 

library(tidyverse)
library(janitor)
library(ggplot2)
library(lubridate)

rm(list = ls(all.names = TRUE))

setwd("/Users/clh89/MATLAB/Projects/Kenya_treatment/KECEA/")  # *******SET ME***********
numSces = 0 # number of scenarios to run through

# Translating Matlab indexes for compartments into R factors
ageCateg = data.frame("index" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17), 
                     "ageCateg" = c("age.0.4", "age.5.9", "age.10.14", "age.15.19", "age.20.24", "age.25.29", "age.30.34", "age.35.39", 
                                      "age.40.44", "age.45.49", "age.50.54", "age.55.59", "age.60.64", "age.65.69", "age.70.74", "age.75.79", "all.ages"))
hpvCateg = data.frame("index" = seq(1,7,1), 
                        "hpvCateg" = c("hpv.susceptible", "hpv.infected", "cin1", "cin2", "cin3", "cc", "hpv.immune"))
dxCateg = data.frame("index" = seq(1,10,1), 
                     "dxCateg" = c("undiagnosed", "undiagnosed", "undiagnosed", "dx, untreated", "dx, untreated", "dx, untreated", "dx, treated", "dx, treated", "dx, treated", "hysterectomy"))
ccStageCateg = data.frame("index" = seq(1,10,1), 
                          "stageCateg" = c("local", "regional", "distant", "local", "regional", "distant", "local", "regional", "distant", "hysterectomy"))
deathsCateg = data.frame("index" = seq(1, 3, 1), 
                         "deathsCateg" = c("treated cc", "untreated cc", "all cause"))
screenTreatCateg = data.frame("index" = seq(1,3,1), 
                              "screenTreatCateg" = c("nScreened", "nColpo", "nCinTreat"))
dxCCTreatCateg = data.frame("treat" = seq(1,3,1), 
                                    "dxType" = c("dx, treated", "dx, untreated", "hysterectomy"))
screenSympCCTreatCateg = data.frame("index" = c(1, 2), 
                                    "screenSymp" = c("screening", "symptoms"))

# loop through each of the scenarios

for (sceNum in seq(0, numSces, 1)) {

    deaths = read.csv(paste0("deaths_S", sceNum, ".csv"))
    screenTreat = read.csv(paste0("screenTreat_S", sceNum, ".csv"))
    hpvHealthState = read.csv(paste0("hpvHealthState_S", sceNum, ".csv"))
    totalPerAge = read.csv(paste0("totalPerAge_S", sceNum, ".csv"))
    vax = read.csv(paste0("vax_S", sceNum, ".csv"))
    ccHealthState = read.csv(paste0("ccHealthState_S", sceNum, ".csv"))
    newCC = read.csv(paste0("newCC_S", sceNum, ".csv"))
    screenSympCCTreat = read.csv(paste0("screenSympCCTreat_S", sceNum, ".csv"))
    
    # DEATHS 
    
    deaths_clean <- deaths %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      left_join(., (totalPerAge %>% rename(N=count)), by=c("sceNum", "paramNum", "year", "age")) %>% 
      mutate(deathCateg = case_when(index == 1 ~ "ccDeath_treat", 
                               index == 2 ~ "ccDeath_untreat", 
                               index == 3 ~ "allDeath"), 
             sceNum = sceNum - 1) %>% # the scenario actually starts at 0, but matlab indices start at 1
      select(sceNum, paramNum, year, ageCateg, N, deathCateg, count) %>% unique %>% 
      spread(., deathCateg, count) %>% 
      mutate(ccDeath = ccDeath_treat + ccDeath_untreat) %>% 
      select(-c(ccDeath_treat, ccDeath_untreat))
    
    # SCREENING AND TREATMENT
    
    screenTreat_clean <- screenTreat %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      left_join(., screenTreatCateg, by=c("index")) %>% 
      mutate(sceNum = sceNum - 1) %>% 
      select(sceNum, paramNum, year, ageCateg, screenTreatCateg, count) %>% unique %>%
      spread(., screenTreatCateg, count)
    
    # HPV HEALTH STATES
    
    hpvHealthState_clean <- hpvHealthState %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      left_join(., hpvCateg, by=c("healthState" = "index")) %>% 
      mutate(sceNum = sceNum - 1) %>% 
      select(sceNum, paramNum, year, ageCateg, hpvCateg, count) %>% unique %>% 
      spread(., hpvCateg, count)
    
    # CERVICAL CANCER HEALTH STATES 
    
    ccHealthState_clean <- ccHealthState %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      left_join(., ccStageCateg, by=c("endpoint" = "index")) %>% 
      left_join(., dxCateg, by=c("endpoint" = "index")) %>% 
      mutate(sceNum = sceNum - 1, 
             stageCateg_ab = case_when(stageCateg == "local" ~ "l", 
                                       stageCateg == "regional" ~ "r", 
                                       stageCateg == "distant" ~ "d"), 
             dxCateg_ab = case_when(dxCateg == "undiagnosed" ~ "udx", 
                                    dxCateg == "dx, treated" ~ "dx", 
                                    dxCateg == "dx, untreated" ~ "ut"), 
             stageDxCateg = case_when(stageCateg == "hysterectomy" ~ "cc_prev_hyst", 
                                      TRUE ~ paste0("cc_", stageCateg_ab, "_prev_", dxCateg_ab))) %>% 
      select(sceNum, paramNum, year, ageCateg, stageDxCateg, count) %>% unique %>% 
      spread(., stageDxCateg, count) 
    
    # CERVICAL CANCER INCIDENCE (DIAGNOSED CASES) STRATIFIED BY TREATMENT TYPE
    
    screenSympCCTreat_clean <- screenSympCCTreat %>% 
          left_join(., ageCateg, by=c("age" = "index")) %>% 
          left_join(., ccStageCateg, by=c("endpoint" = "index")) %>% 
          left_join(., dxCCTreatCateg, by=c("treat")) %>% 
          left_join(., screenSympCCTreatCateg, by=c("index")) %>% 
          mutate(sceNum = sceNum - 1, 
                 stageCateg_ab = case_when(stageCateg == "local" ~ "l", 
                                           stageCateg == "regional" ~ "r", 
                                           stageCateg == "distant" ~ "d"), 
                 dxCateg_ab = case_when(dxType == "dx, treated" ~ "dx", 
                                        dxType == "dx, untreated" ~ "ut"),  
                 stageDxCateg = case_when(dxType == "hysterectomy" ~ "cc_inc_hyst", 
                                          TRUE ~ paste0("cc_", stageCateg_ab, "_inc_", dxCateg_ab))) %>% 
          group_by(sceNum, paramNum, year, ageCateg, stageDxCateg) %>% 
                  mutate(count2 = sum(count, na.rm=TRUE)) %>% ungroup %>% 
          select(paramNum, sceNum, year, ageCateg, stageDxCateg, count2) %>% unique %>% 
          spread(., stageDxCateg, count2) 
    
    # NEW CERVICAL CANCER CASES 
    
    newCC_clean <- newCC %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      mutate(sceNum = sceNum - 1) %>% 
      select(sceNum, paramNum, year, ageCateg, count) %>% unique %>% 
      rename(newCC = count)
    
    # VACCINATIONS
    vax_clean <- vax %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      mutate(sceNum = sceNum - 1, 
             vaxType = case_when(vaxType == 1 ~ "N_vax_school", 
                                 vaxType == 2 ~ "N_vax_cu")) %>% 
      select(sceNum, paramNum, year, ageCateg, vaxType, count) %>% unique %>% 
      spread(., vaxType, count) 
    
    # Combine all the dfs
    combined <- deaths_clean %>% 
      left_join(., screenTreat_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      left_join(., hpvHealthState_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      left_join(., ccHealthState_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      left_join(., screenSympCCTreat_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      left_join(., vax_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      left_join(., newCC_clean, by=c("sceNum", "paramNum", "year", "ageCateg")) %>% 
      select(sceNum, paramNum, year, ageCateg, N, hpv.susceptible, hpv.immune, hpv.infected, cin1, cin2, cin3,
              cc_l_prev_udx, cc_l_prev_ut, cc_l_prev_dx, 
              cc_r_prev_udx, cc_r_prev_ut, cc_r_prev_dx,
              cc_d_prev_udx, cc_d_prev_ut, cc_d_prev_dx, 
              cc_prev_hyst, 
              cc_l_inc_ut, cc_l_inc_dx, 
              cc_r_inc_ut, cc_r_inc_dx, 
              cc_d_inc_ut, cc_d_inc_dx,
              cc_inc_hyst, 
              newCC, nScreened, nCinTreat, nColpo, 
              allDeath, ccDeath, N_vax_school, N_vax_cu) %>% unique
    
  sceDf <- combined %>% 
    filter(sceNum == sceNum)
  
  ifelse(!dir.exists(paste0(getwd(), "/Outputs")), dir.create(paste0(getwd(), "/Outputs")), FALSE)
  
  write.csv(sceDf, paste0(getwd(), "/Outputs/modelResultsForCea_S", sceNum, ".csv"), row.names = FALSE)
}