# title: cleanMatlabOutputs.R
# date: 6/2/2023
# author: Christine L Hathaway
# description: FOR KENYA MODEL. Takes the CSV files spit out by Matlab (vaxCEA_multSims_mainFunction.m) and puts it into a format 
# that can be used for the economic part of the analysis. 


####### TODO
# - Make correction for the first year of vax (need to decide how/if we want to address this)
# - Make a correction for 2 vs 1 dose; will depend on the scenario (I added, need to debug)
# - Make a correction for the 0.7/0.9 coverage adjustment. You have to multiply by 0.9/0.7 to get the actual coverage value. (I added, need to debug)

library(tidyverse)
library(janitor)
library(ggplot2)
library(lubridate)

rm(list = ls(all.names = TRUE))

setwd("/Users/clh89/MATLAB/Projects/Kenya_treatment/KECEA/")  # *******SET ME***********
# numSces = 0 # number of scenarios to run through
sces = c(0, 1, 2, 3, 4, 5, 6, 7, 8)

# Translating Matlab indexes for compartments into R factors
ageCateg = data.frame("index" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17), 
                     "ageCateg" = c("age.0.4", "age.5.9", "age.10.14", "age.15.19", "age.20.24", "age.25.29", "age.30.34", "age.35.39", 
                                      "age.40.44", "age.45.49", "age.50.54", "age.55.59", "age.60.64", "age.65.69", "age.70.74", "age.75.79", "all.ages")) %>% 
                  mutate(ageCateg = factor(ageCateg, levels=c("age.0.4", "age.5.9", "age.10.14", "age.15.19", "age.20.24", "age.25.29", "age.30.34", "age.35.39", 
                                                              "age.40.44", "age.45.49", "age.50.54", "age.55.59", "age.60.64", "age.65.69", "age.70.74", "age.75.79", "all.ages")))
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
sceDose = data.frame("index" = seq(0, 16, 1), 
                     "dose" = c(2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

# loop through each of the scenarios

for (sceNum in sces) {

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
    # TODO: decide whether or not we want to address the issue with the vax burn in period
    vax_clean <- vax %>% 
      left_join(., ageCateg, by=c("age" = "index")) %>% 
      mutate(sceNum = sceNum - 1, 
             vaxType = case_when(vaxType == 1 ~ "N_vax_school", 
                                 vaxType == 2 ~ "N_vax_cu")) %>% 
      select(sceNum, paramNum, year, ageCateg, vaxType, count) %>% unique %>% 
      spread(., vaxType, count) %>% 
      # for the vaccines in year 2023, the first timepoint is a burn in period, and is not reflective of the actual number of vaccines that would have been
      # administered. so just for 2023.00, i take the num vaccines from the timepoint after. 
      # don't need to do this for the catchup vaxxes because it is a muli-age-cohort vaccination
          
      # 6/13/23: I am leaning towards not doing this anymore. In the grand scheme of things, if there is a spike in 2023, these 5-8 year olds will be vaccinated at some point
      # we are just artificially bringing it forward in time to get around the 5-year age group issue. 
      # group_by(sceNum, paramNum, ageCateg) %>% 
      #       arrange(sceNum, paramNum, ageCateg, year) %>% 
      #       mutate(N_vax_school_2 = case_when(year == 2023.167 & ageCateg == "all.ages" ~ lead(N_vax_school), 
      #                                       TRUE ~ N_vax_school)) %>% ungroup %>% 
      #     filter(paramNum == 18, ageCateg == "all.ages") %>% View()
      # correct for coverage adjustment for quadrivalent (0.7/0.9)
      mutate(N_vax_school = N_vax_school * (0.9/0.7), 
             N_vax_cu = N_vax_cu * (0.9/0.7)) %>% 
      # times 2 for 2-dose scenarios
      left_join(sceDose, by=c("sceNum"="index")) %>%
      mutate(N_vax_school = case_when(dose == 2 ~ N_vax_school * 2, 
                                      dose == 1 & year>= 2019 & year <2023 ~ N_vax_school * 2, 
                                      TRUE ~ N_vax_school), 
             N_vax_cu = case_when(dose == 2 ~ N_vax_cu * 2, 
                                  dose == 1 & year>= 2019 & year <2023 ~ N_vax_cu * 2, 
                                  TRUE ~ N_vax_cu)) %>%  
      select(sceNum, paramNum, year, ageCateg, N_vax_cu, N_vax_school) 
    
    # two years that matter
    
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
  
  print(paste0("Scenario ", sceNum, " complete."))
}

######### Fixing issue where the scenario numbers in the CSV don't actually match with the actual scenario number #########

scenarios <- seq(0,8,1)

for (sce in scenarios) {
      
      read <- read.csv(paste0("/Users/clh89/MATLAB/Projects/Kenya_treatment/KECEA/Outputs/modelResultsForCea_S", sce, ".csv"))
      
      read <- read %>% mutate(sceNum = sce)
      
      write.csv(read, paste0(getwd(), "/Outputs/modelResultsForCea_S", sce, ".csv"), row.names = FALSE)
}