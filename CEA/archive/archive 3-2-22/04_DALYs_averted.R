# DALYs
# October 28, 2020


## REVISIONS REQUIRED:

# 1. INCLUDE ALL-CAUSE MORTALITY!!!
# 2. NEED TO DISCOUNT!!! 
# 3. THESE CALCULATIONS MUST BE REVISED TO REFLECT FOX-RUSHBY calculations !!
# 4. Note that this code extends the calculation for half-years to the full year.

# Note that we choose NOT to age-weight (as in the original DALY calculations by Murray /Fox-Rushby et al) because
# we want to weight all lives equally

########################################################################################################################

# Set up disability weight parameters

wt_hiv_neg <- 0
wt_cd4_350 <- 0.078
wt_cd4_250_350 <- 0.274
wt_cd4_200 <- 0.582
wt_art <- 0.078
wt_dead <- 1

cd4 <- c("CD4_350plus_noART",  "CD4_200-350_noART", "CD4_below200_noART","onART")
wt <- c(wt_cd4_350, wt_cd4_250_350, wt_cd4_200, wt_art)

wtDF <- as.data.frame( cbind(cd4, wt)) 

wt_shiftDF <- wtDF %>% filter(cd4!="onART") %>% 
  mutate(wt_shift = as.character( as.numeric(wt) - lag(as.numeric(wt), default = 0))) %>% 
  filter(cd4!="CD4_350plus_noART") %>% 
  mutate(cd4= substr(cd4, 1, nchar(cd4)-6))

wt_all <- c("HIV_neg", "CD4_350plus_noART",  "CD4_200-350_noART", "CD4_below200_noART","onART")
wt_CD4_all <- c(wt_cd4_350, wt_cd4_250_350, wt_cd4_200, wt_art)

###########################################################################################################################

# Set up age dataframe

timeHorizon <- horizon_year   # (or add 1? )We count health effects until the END of 2060
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


#################################################################################################################

for (i in 1:3) {   # SCENARIO

  #########################
  # INCIDENCE, by age bin #  
  #########################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_incidence_combined_",a,".csv"), header=F) %>%  
      
      # Rename columns (29 columns: year, mean, min, max, 25 sets)
      setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
      
      # Calculate Life-Years Left (truncated at the year of the Time Horizon)
      mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
             yearStop = year + yearsToEnd) 
    
    # Important to discount health effects by year, which requires separating out the year of the health effect, rather than the year of the incident
    # To address this, generate a new DF for each year of incidents, then multiply by the weight, discount, and aggregate
    
    for (y in unique(x$year)) {
      
      yearStop <- x[x$year==y,"yearStop"]
      yearStop_rounded <- floor(x[x$year==y,"yearStop"]) # rounds down so that years are whole numbers 
      
      # What years are excluded from this DF
      if ( y > 2020 & yearStop_rounded < timeHorizon)  {otherYears <- data.frame( year= c(2020:(y-1), (yearStop_rounded+1):timeHorizon))} 
      if ( y == 2020  & yearStop_rounded < timeHorizon) {otherYears <- data.frame( year = (yearStop_rounded+1):timeHorizon) }
      if ( y > 2020 & yearStop_rounded == timeHorizon) {otherYears <- data.frame( year = c(2020:(y-1))) }
      if ( y == 2020 & yearStop_rounded == timeHorizon) {otherYears <- NULL }
      
      if (timeHorizon  - y > yearStop_rounded + y )  {otherYears <- data.frame( year = (yearStop_rounded+1):timeHorizon) }
      
      v <- x[x$year==y, -c(1,30,31)] # keep only our 28 parameter values
      
      # Multiply incidence by the  Weight for cases with CD4+ >350, bind with rows of 0's for previous years
      newDF <- data.frame(year = y:yearStop_rounded) %>%
        cbind(v * wt_cd4_350) %>%   
        plyr::rbind.fill(otherYears) %>%  # fills empty rows with NA
        replace(is.na(.), 0) %>%   # replace NAs with 0's
        arrange(year) %>% select(-year) %>% 
        
        # DISCOUNT
        discount(., discount_rate = discount_rate) %>%  
        
        addYearCol(., horizon_year = horizon_year) %>% 
        recalcFuns(.) # recalculate mean, min, max) 
        
      assign(paste0("temp",y), newDF)                    
     
    }
    
    t <- mget(ls(pattern="temp"))
    TEMP <- Reduce(`+`, t)   %>%  # SUM Across ALL DFs
      mutate(year=2020:2060) %>%  select(year, everything()) 
      
    assign(paste0("incScen",i,"_",a), TEMP)
    
    rm(list=ls(pattern="temp"))
    rm(v, t, TEMP, newDF)
  }
  
  # Sum across ages, Save, Housekeeping 
  
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
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year)) %>%  
        
        # Calculate Life-Years Left (truncated at the year of the Time Horizon)
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
             yearStop = year + yearsToEnd) 
      
      # Important to discount health effects by year, which requires separating out the year of the health effect, rather than the year of the incident
      # To address this, generate a new DF for each year of incidents, then multiply by the weight, discount, and aggregate
      
      for (y in unique( x$year )) {
        
        yearStop <- x[x$year==y,"yearStop"]
        yearStop_rounded <- floor(x[x$year==y,"yearStop"]) # rounds down so that years are whole numbers 
        
        # What years are excluded from this DF
        if ( y > 2020 & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year= c(2020:(y-1), (yearStop_rounded+1):timeHorizon))
          } 
        if ( y == 2020  & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year = (yearStop_rounded+1):timeHorizon) 
          }
        if ( y > 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- data.frame( year = c(2020:(y-1))) 
          }
        if ( y == 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- NULL 
          }
      
        v <- x[ x$year==y, -c(1,30,31)] # keep only our 28 parameter values
        
        # Multiply incidence by the  corresponding SHIFT in utility weight, bind with rows of 0's for previous years
        newDF <- data.frame(year = y:yearStop_rounded) %>%
          cbind(v * wt_shift) %>%   
          plyr::rbind.fill(otherYears) %>%  # fills empty rows with NA
          replace(is.na(.), 0) %>%   # replace NAs with 0's
          arrange(year) %>% select(-year) %>% 
          
          # DISCOUNT
          discount(., discount_rate = discount_rate) %>%  
          
          addYearCol(., horizon_year = horizon_year) %>% 
          recalcFuns(.) # recalculate mean, min, max
          
          addYearCol(., horizon_year = horizon_year) 
        
        assign(paste0("temp",y), newDF)                    
        
      }
      
      t <- mget(ls(pattern="temp"))
      TEMP <- Reduce(`+`, t)   %>%  # SUM Across ALL DFs
        mutate(year=2020:2060) %>%  select(year, everything()) 
      
      assign(paste0("cd4ShiftScen",i,"_",c,"_",a), TEMP)
      
      rm(list=ls(pattern="temp"))
      rm(v, t, TEMP, newDF)
    }
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="cd4Shift"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYcd4_Scen",i), DALY)
  
  rm(list=ls(pattern="cd4Shift"))
  rm(x,d, DALY)
  
  
  #################################################
  # HIV-Associated Deaths, by CD4+ count, age bin #    -- NEED TO SWITCH THIS TO FULL HIV-POSITIVE DEATHS, AND ADD BACKGROUND HIV-NEGATIVE DEATHS
  #################################################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in cd4) {  # CD4+ CATEGORY
      
      wt_cd4 <- wtDF[wtDF$cd4==c, "wt"]           # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_mortality_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year))
      
      # Important to discount health effects by year, which requires separating out the year of the health effect, rather than the year of the incident
      # To address this, generate a new DF for each year of incidents, then multiply by the weight, discount, and aggregate
      
      for (y in unique(x$year)) {
        
        yearStop <- x[x$year==y,"yearStop"]
        yearStop_rounded <- floor(x[x$year==y,"yearStop"]) # rounds down so that years are whole numbers 
        
        # What years are excluded from this DF
        if ( y > 2020 & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year= c(2020:(y-1), (yearStop_rounded+1):timeHorizon))
        } 
        if ( y == 2020  & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year = (yearStop_rounded+1):timeHorizon) 
        }
        if ( y > 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- data.frame( year = c(2020:(y-1))) 
        }
        if ( y == 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- NULL 
        }
        
        v <- x[x$year==y, -c(1,30,31)] # keep only our 28 parameter values
        
        # Multiply incidence by the  corresponding SHIFT in utility weight, bind with rows of 0's for previous years
        
        newDF <- data.frame(year = y:yearStop_rounded) %>%
          
          # APPLY UTILITY WEIGHT
          cbind(v * (wt_dead - as.numeric(wt_cd4)) ) %>%      
          plyr::rbind.fill(otherYears) %>%  # fills empty rows with NA
          replace(is.na(.), 0) %>%   # replace NAs with 0's
          arrange(year) %>% select(-year) %>% 
          
          # DISCOUNT
          discount(., discount_rate = discount_rate) %>%  
          
          addYearCol(., horizon_year = horizon_year) %>% 
          recalcFuns(.) # recalculate mean, min, max
        
        assign(paste0("temp",y), newDF)                    
        
      }
      
      t <- mget(ls(pattern="temp"))
      TEMP <- Reduce(`+`, t)   %>%  # SUM Across ALL DFs
        mutate(year=2020:2060) %>%  select(year, everything()) 
      
      assign(paste0("HIVmortScen",i,"_",c,"_",a), TEMP)
      
      rm(list=ls(pattern="temp"))
      rm(v, t, TEMP, newDF)
      
    }
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="HIVmort"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYdeathHIV_Scen",i),DALY)
  
  rm(list=ls(pattern="HIVmort"))
  rm(x, d, DALY)
  
  #################################################
  # BACKGROUND Deaths, by CD4+ count, age bin #   -- NEED TO UPDATE
  #################################################
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in cd4) {  # CD4+ CATEGORY
      
      wt_cd4 <- wtDF[wtDF$cd4==c, "wt"]           # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_mortality_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year))
      
      # Important to discount health effects by year, which requires separating out the year of the health effect, rather than the year of the incident
      # To address this, generate a new DF for each year of incidents, then multiply by the weight, discount, and aggregate
      
      for (y in unique(x$year)) {
        
        yearStop <- x[x$year==y,"yearStop"]
        yearStop_rounded <- floor(x[x$year==y,"yearStop"]) # rounds down so that years are whole numbers 
        
        # What years are excluded from this DF
        if ( y > 2020 & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year= c(2020:(y-1), (yearStop_rounded+1):timeHorizon))
        } 
        if ( y == 2020  & yearStop_rounded < timeHorizon) {
          otherYears <- data.frame( year = (yearStop_rounded+1):timeHorizon) 
        }
        if ( y > 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- data.frame( year = c(2020:(y-1))) 
        }
        if ( y == 2020 & yearStop_rounded == timeHorizon) {
          otherYears <- NULL 
        }
        
        v <- x[x$year==y, -c(1,30,31)] # keep only our 28 parameter values
        
        # Multiply incidence by the  corresponding SHIFT in utility weight, bind with rows of 0's for previous years
        
        newDF <- data.frame(year = y:yearStop_rounded) %>%
          # APPLY UTILITY WEIGHT
          cbind(v * (wt_dead - as.numeric(wt_cd4)) ) %>%      
          plyr::rbind.fill(otherYears) %>%  # fills empty rows with NA
          replace(is.na(.), 0) %>%   # replace NAs with 0's
          arrange(year) %>% select(-year) %>% 
          
          # DISCOUNT
          discount(., discount_rate = discount_rate) %>%  
          
          addYearCol(., horizon_year = horizon_year) %>% 
          recalcFuns(.) # recalculate mean, min, max
        
        assign(paste0("temp",y), newDF)                    
        
      }
      
      t <- mget(ls(pattern="temp"))
      TEMP <- Reduce(`+`, t)   %>%  # SUM Across ALL DFs
        addYearCol(., horizon_year = horizon_year)
      
      assign(paste0("HIVmortScen",i,"_",c,"_",a), TEMP)
      
      rm(list=ls(pattern="temp"))
      rm(v, t, TEMP, newDF)
      
    }
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="HIVmort"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    addYearCol(., horizon_year = horizon_year) 
  
  assign(paste0("DALYdeathHIV_Scen",i),DALY)
  
  rm(list=ls(pattern="HIVmort"))
  rm(x, d, DALY)
  
}

###################################################################################################################

# Annual

daly_scen1 <- (DALYcases_Scen1[,-1] + DALYcd4_Scen1[,-1] + DALYdeathHIV_Scen1[,-1] + DALYdeathBKGD_Scen1[,-1]) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_scen2 <- (DALYcases_Scen2[,-1] + DALYcd4_Scen2[,-1] + DALYdeathHIV_Scen2[,-1] + DALYdeathBKGD_Scen2[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_scen3 <- (DALYcases_Scen3[,-1] + DALYcd4_Scen3[,-1] + DALYdeathHIV_Scen3[,-1] + DALYdeathBKGD_Scen3[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative 

for (x in 1:3) {
  
  daly <- get(paste0("daly_scen",x))[,-1] %>% 
    transmute_at(1:28, ~cumsum(.)) %>% 
    addYearCol(., horizon_year = horizon_year) 
  
  assign(paste0("daly_cum_scen",x),daly)
  
}

##################################################################################################################

# Annual raw dalys averted 

daly_averted_raw_scen2 <- (daly_scen1[,-1] - daly_scen2[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_averted_raw_scen3 <- (daly_scen1[,-1] - daly_scen3[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Annual percent dalys averted

daly_averted_pct_scen2 <- (((daly_scen1[,-1] - daly_scen2[,-1] )/ daly_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_averted_pct_scen3 <- (((daly_scen1[,-1] - daly_scen3[,-1] )/ daly_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative raw dalys averted 

daly_cum_averted_raw_scen2 <- ( daly_cum_scen1[,-1] - daly_cum_scen2[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_cum_averted_raw_scen3 <- ( daly_cum_scen1[,-1] - daly_cum_scen3[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative percent dalys averted

daly_cum_averted_pct_scen2 <- (((daly_cum_scen1[,-1] - daly_cum_scen2[,-1] )/ daly_cum_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_cum_averted_pct_scen3 <- (((daly_cum_scen1[,-1] - daly_cum_scen3[,-1] )/ daly_cum_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

###################################################################################################################

# Export CSVs

csv_list <- list("daly_scen1",
                 "daly_scen2",
                 "daly_scen3",
                 "daly_cum_scen1",
                 "daly_cum_scen2",
                 "daly_cum_scen3",
                 "daly_averted_raw_scen2",
                 "daly_averted_raw_scen3",
                 "daly_averted_pct_scen2",
                 "daly_averted_pct_scen3",
                 "daly_cum_averted_raw_scen2",
                 "daly_cum_averted_raw_scen3",
                 "daly_cum_averted_pct_scen2",
                 "daly_cum_averted_pct_scen3")

lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"effects/daly/",x,".csv"), row.names = F))


####################################################################################################################

# Remove excess DFs, source for next script

rm(list=ls()[! ls() %in% c("main_path","cea_path","helper_path")])
source(paste0(helper_path, "helper.R"))

