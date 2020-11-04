# DALYs
# October 28, 2020

library(tidyverse)

########################################################################################################################

# Set up disability weight parameters

wt_hiv_neg <- 0
wt_cd4_350 <- 0.078
wt_cd4_250_350 <- 0.274
wt_cd4_200 <- 0.582
wt_art <- 0.078
wt_dead <- 1

cd4 <- c("CD4_350plus_noART",  "CD4_200-350_noART", "CD4_below200_noART","onART")
wt <- c (wt_cd4_350, wt_cd4_250_350, wt_cd4_200, wt_art)

wtDF <- as.data.frame( cbind(cd4, wt)) 

wt_shiftDF <- wtDF %>% filter(cd4!="onART") %>% 
  mutate(wt_shift = as.character( as.numeric(wt) - lag(as.numeric(wt), default = 0))) %>% 
  filter(cd4!="CD4_350plus_noART") %>% 
  mutate(cd4= substr(cd4, 1, nchar(cd4)-6))

###########################################################################################################################

# Set up age dataframe

timeHorizon <- 2060 + 1   # We count until the END of 2060
lifeExp <- 80

ages = c("15-19",
        "20-24",
        "25-29",
        "30-34",
        "35-39",
        "40-44",
        "45-49",
        "50-54",
        "55-59",
        "60-64",
        "65-69",
        "70-74",
        "75-79")

ageDF <- data.frame( ages = ages,
                     ages.midpt = c(seq(17.5,77.5,5)) )%>% 
  mutate(yearsExp = lifeExp-ages.midpt)

##########################################################################################################################

# CALCULATE DALYs

for (i in 1:2) {   # SCENARIO

  #########################
  # INCIDENCE, by age bin #  -- CONFIRM WITH CARA CORRECT FILES !!!!
  #########################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_incidence_combined_",a,".csv"), header=F) %>%  
      
      # Housekeeping 
      setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
      
      # Calculate Life-Years Left
      mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year)) %>%  
      
      # Multiply # New Cases by Life Years Left  * Weight for cases with CD4+ >350
      transmute_at(2:29, ~(.) * yearsToEnd * (wt_cd4_350) )     # MAKE SURE OK TO ASSUME CD4+ > 350 IS OK
    
    assign(paste0("incScen",i,"_",a), x)
  }
  
  # Sum, Save, Housekeeping 
  
  d <- mget(ls(pattern="inc"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYcases_Scen",i), DALY)
  
  rm(list=ls(pattern="inc"))
  rm(x, d, DALY)
  
  ##############################
  # CD4+ count shifts, age bin #
  ##############################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in unique(wt_shiftDF$cd4)) {  # CD4+ CATEGORY
      
      wt_shift <- wt_shiftDF[wt_shiftDF$cd4==c, "wt_shift"]           # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_CD4_incidence_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Housekeeping 
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year)) %>%  
        
        # Multiply # Cases by Life Years Left  * (DIFFERENCE in Disability Weight for death and weight for the CD4+ count category at death)
        transmute_at(2:29, ~(.) * yearsToEnd * (as.numeric(wt_shift)) )  
      
      assign(paste0("cd4ShiftScen",i,"_",c,"_",a), x)
    }
  }
  
  # Sum, Save, Housekeeping 
  
  d <- mget(ls(pattern="cd4Shift"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYcd4_Scen",i), DALY)
  
  rm(list=ls(pattern="cd4Shift"))
  rm(x,d, DALY)
  
  
  ##################################
  # Deaths, by CD4+ count, age bin #
  ##################################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in cd4) {  # CD4+ CATEGORY
      
      wt_cd4 <- wtDF[wtDF$cd4==c, "wt"]           # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_mortality_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Housekeeping 
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year)) %>%  
        
        # Multiply # Deaths by Life Years Left  * (DIFFERENCE in Disability Weight for death and weight for the CD4+ count category at death)
        transmute_at(2:29, ~(.) * yearsToEnd * (wt_dead - as.numeric(wt_cd4) ) )  
      
      assign(paste0("mortScen",i,"_",c,"_",a), x)
      
    }
    
  }
  
  # Sum, Save, Housekeeping 
  
  d <- mget(ls(pattern="mort"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYdeath_Scen",i),DALY)
  
  rm(list=ls(pattern="mort"))
  rm(x, d, DALY)
  
}







