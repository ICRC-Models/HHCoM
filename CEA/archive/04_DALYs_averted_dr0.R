# DALYs
# January 13, 2021

# Note that we choose NOT to age-weight (as in the original DALY calculations by Murray /Fox-Rushby et al) because
# we want to weight all lives equally

# 1. ART initiation and discontinutation
# 2. Update discounting to year of health effect
# 3. Generate graphs - DIFFERENCE in DALYs between the Scenarios, and disaggregation

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

cd4_2 <- c("CD4_350plus",  "CD4_200-350", "CD4_below200")
wt_2 <- c(wt_cd4_350, wt_cd4_250_350, wt_cd4_200)
wt_ART <- as.data.frame(cbind(cd4_2, wt_2)) 

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

#########################
# INCIDENCE, by age bin #  
#########################

for (i in 1:3) {   # SCENARIO
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_HIV_incidence_combined_",a,".csv"), header=F) %>%  
      
      # Rename columns (29 columns: year, mean, min, max, 25 sets)
      setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
      
      # Calculate Life-Years Left (truncated at the year of the Time Horizon) 
      mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year))  
      
      # APPLY UTILIY WEIGHT for HIV * years at that weight
     newDF <- ( x[ ,2:29] * x$yearsToEnd * wt_cd4_350 )  
      
      assign(paste0("inc",a), newDF)                    
      
    }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="inc"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year = 2020:2060) %>% select(year, everything()) 
  
  assign(paste0("DALYcases_Scen",i), DALY)
  
  rm(list=ls(pattern="inc"))
  rm(x, d, DALY)
  
}

##############################
# CD4+ count shifts, age bin #
##############################

for (i in 1:3) {   # SCENARIO  
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in unique(wt_shiftDF$cd4)) {  # 2 CD4+ CATEGORYs
      
      wt_shift <- wt_shiftDF[wt_shiftDF$cd4==c, "wt_shift"]   # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_CD4_incidence_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left (truncated at the year of the Time Horizon)
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
               yearStop = year + yearsToEnd) 
      
        # Multiply incidence by the  corresponding SHIFT in utility weight
        newDF <- ( x[ ,2:29] * x$yearsToEnd * as.numeric(wt_shift) )  
        
        assign(paste0("temp",c), newDF)                    
        
      }
    
    # Sum across weights, Save, Housekeeping 
    
    w <- mget(ls(pattern="temp"))
    
    TEMP <- Reduce(`+`, w)   %>%  # SUM Across ALL DFs
      mutate(year = 2020:2060) %>% select(year, everything()) 
    
    assign(paste0("cd4Shift_Scen",a), TEMP)
    
   rm(list=ls(pattern="temp"))
   rm(x,w, TEMP)
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="cd4Shift"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year = 2020:2060) %>% select(year, everything()) 
  
  assign(paste0("DALYcd4_Scen",i), DALY)
  
 rm(list=ls(pattern="cd4Shift"))
 rm(d, DALY)
  
}    

#########################################
# ART initiation by CD4+ count, age bin #
#########################################

# Remember you've deleted the discounting here
# THis is given weird results?? 

for (i in 1:3) {   # SCENARIO  
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in unique(wt_shiftDF$cd4)) {  # CD4+ CATEGORY
      
      wt_shift <- wt_shiftDF[wt_shiftDF$cd4==c, "wt_shift"] # Weight ADDED by starting on ART
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_ART_incidence_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left (truncated at the year of the Time Horizon)
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
               yearStop = year + yearsToEnd) 
      
      # Multiply incidence by the  corresponding SHIFT in utility weight
      newDF <- ( x[ ,2:29] * x$yearsToEnd * as.numeric(wt_shift) * (-1) )  # Multiply by (-1) because adding it
      
      assign(paste0("temp",c), newDF)                    
      
    }
    
    # Sum across weights, Save, Housekeeping 
    
    w <- mget(ls(pattern="temp"))
    
    TEMP <- Reduce(`+`, w)   %>%  # SUM Across ALL DFs
      mutate(year = 2020:2060) %>% select(year, everything()) 
    
    assign(paste0("artShift_Scen", a), TEMP)
    
    rm(list=ls(pattern="temp"))
    rm(x,w, TEMP)
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="artShift"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year = 2020:2060) %>% select(year, everything()) 
  
  assign(paste0("DALYart_Scen",i), DALY)
  
  rm(list=ls(pattern="artShift"))
  rm(d, DALY)
  
}    

###############################################
# ART discontinuation by CD4+ count, age bin #
###############################################

# NOT DISCOUNTED

for (i in 1:3) {   # SCENARIO  
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in unique(wt_shiftDF$cd4)) {  # CD4+ CATEGORY
      
      wt_shift <- wt_shiftDF[wt_shiftDF$cd4==c, "wt_shift"]   # Weight LOST by stopping ART
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_ART_discont_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left (truncated at the year of the Time Horizon)
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
               yearStop = year + yearsToEnd) 
      
      # Multiply incidence by the  corresponding SHIFT in utility weight
      newDF <- ( x[ ,2:29] * x$yearsToEnd * as.numeric(wt_shift) )  
      
      assign(paste0("temp",c), newDF)                    
      
    }
    
    # Sum across weights, Save, Housekeeping 
    
    w <- mget(ls(pattern="temp"))
    
    TEMP <- Reduce(`+`, w)   %>%  # SUM Across ALL DFs
      mutate(year = 2020:2060) %>% select(year, everything()) 
    
    assign(paste0("artDiscShift_Scen", a), TEMP)
    
    rm(list=ls(pattern="temp"))
    rm(x,w, TEMP)
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="artDiscShift"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year = 2020:2060) %>% select(year, everything()) 
  
  assign(paste0("DALYartDISC_Scen",i), DALY)
  
  rm(list=ls(pattern="artDiscShift"))
  rm(d, DALY)
  
}    


##############################################################
# All-Cause HIV-Positive Deaths, by ART, CD4+ count, age bin #    
##############################################################

for (i in 1:3) {   # SCENARIO  
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    for (c in cd4) {  # CD4+ CATEGORY
      
      wt_cd4 <- wtDF[wtDF$cd4==c, "wt"]           # Pull the corresponding utility weight for that CD4+ cat
      
      x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_allCause_mortality_combined_",a,"_",c,".csv"), header=F) %>%  
        
        # Rename columns (29 columns: year, mean, min, max, 25 sets)
        setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
        
        # Calculate Life-Years Left (truncated at the year of the Time Horizon)
        mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year),
               yearStop = year + yearsToEnd) 
    
        # Multiply  mortality by the  corresponding SHIFT in UTILITY WEIGHT
        
        newDF <- ( x[ ,2:29] * x$yearsToEnd * (wt_dead - as.numeric(wt_cd4))) 
        
        assign(paste0("temp",c), newDF)                    
        
      }
    
    # Sum across weights, Save, Housekeeping 
    
    w <- mget(ls(pattern="temp"))
    
    TEMP <- Reduce(`+`, w)   %>%  # SUM Across ALL DFs
      mutate(year=2020:2060) %>%  select(year, everything()) 
    
    assign(paste0("HIVmort_Scen", a), TEMP)
    
    rm(list=ls(pattern="temp"))
    rm(x, w, TEMP)
    
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="HIVmort"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year=2020:2060) %>%  select(year, everything()) 
  
  assign(paste0("DALYdeathHIV_Scen",i),DALY)
  
  rm(list=ls(pattern="HIVmort"))
  rm(d, DALY)
  
}  


##############################################
# All-cause HIV-negative deaths, by age bin #    
#############################################

for (i in 1:3) {   # SCENARIO
  
  for (a in ages) { # AGE BIN
    
    yearsExp <- ageDF[ageDF$ages==a,"yearsExp"]   # Pull the corresponding years of life remaining for that age
    
    x <- read.csv(paste0(main_path,"Scenario",i,"/Raw_allCause_mortality_combined_",a,"_HIV_neg.csv"), header=F) %>%  
      
      # Rename columns (29 columns: year, mean, min, max, 25 sets)
      setNames(paste0("s",-3:25)) %>% dplyr::rename(year=1, mean=2, min=3, max=4) %>% filter(year>=2020) %>%  
      
      # Calculate Life-Years Left (truncated at the year of the Time Horizon)
      mutate(yearsToEnd = ifelse(timeHorizon-year > yearsExp, yearsExp, timeHorizon-year)) 
      
      # Multiply by full utility weight (1)
         newDF <- ( x[ ,2:29] * x$yearsToEnd * (wt_dead)) %>% 
        
        mutate(year = 2020:2060) %>% select(year, everything()) %>% 
        recalcFuns(.) # recalculate mean, min, max) 
      
      assign(paste0("mort",a), newDF)                    
      
  }
  
  # Sum across ages, Save, Housekeeping 
  
  d <- mget(ls(pattern="mort"))
  
  DALY <- Reduce(`+`, d)   %>%  # SUM Across ALL DFs
    mutate(year = 2020:2060) %>% select(year, everything()) 
  
  assign(paste0("DALYdeathNoHIV_Scen",i), DALY)
  
  rm(list=ls(pattern="mort"))
  rm(x, d, DALY)
  
}
  

###################################################################################################################

# Annual

daly_scen1 <- (DALYcases_Scen1[,-1] + DALYcd4_Scen1[,-1] + DALYart_Scen1[,-1] + DALYartDISC_Scen1[,-1] +
                 DALYdeathHIV_Scen1[,-1] + DALYdeathNoHIV_Scen1[,-1]) %>% 
  # DISCOUNT
 # discount(., discount_rate = discount_rate) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_scen2 <- (DALYcases_Scen2[,-1] + DALYcd4_Scen2[,-1] + DALYart_Scen2[,-1] + DALYartDISC_Scen2[,-1] +
                 DALYdeathHIV_Scen2[,-1] + DALYdeathNoHIV_Scen2[,-1] ) %>% 
  # DISCOUNT
 # discount(., discount_rate = discount_rate) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

daly_scen3 <- (DALYcases_Scen3[,-1] + DALYcd4_Scen3[,-1] + DALYart_Scen3[,-1] + DALYartDISC_Scen3[,-1] +
                 DALYdeathHIV_Scen3[,-1] + DALYdeathNoHIV_Scen3[,-1] ) %>% 
  # DISCOUNT
 # discount(., discount_rate = discount_rate) %>% 
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

lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"results/dr0/daly/",x,".csv"), row.names = F))


####################################################################################################################

# Remove excess DFs, source for next script

rm(list=ls()[! ls() %in% c("main_path","cea_path","helper_path")])
source(paste0(helper_path, "helper.R"))

