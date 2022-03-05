# Calculate QALYs
# MSahu
# September 9, 2020

# QALY Calculation: Years of Life * Utility Value
# This script calculates average annual QALYs gained between interventions


# Note that health benefits for INCIDENT cases are not included. - ie. 
# assumption is that health benefits show up the next year.


##############################################################################


# Set up utility weight parameters

wt_hiv_neg <- 1
wt_cd4_350 <- 0.94
wt_cd4_250_350 <- 0.82
wt_cd4_200 <- 0.7
wt_art <- 0.94
wt_dead <- 0

# Current strategy - calculate on an annual basis ; from the START of the year and don't include incident cases.
#1. Number of prevalent cases * utility weight (by CD4+ count) 
#2. Number of HIV negative cases * utility weight  [can be calculated as 1 - prevalence]
#3. Sum utilities
#4. Discount

# Then Ask Paul Revill - because actually should be cumulative ! (or is there a way to track births?)


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


# MULTIPLY BY UTILITY WEIGHTS 

for (x in 1:3) {
  
  hiv_pos <- (get(paste0("prev_ART_scen",x))[,-1] * wt_art +    
              get(paste0("prev_200_scen",x))[,-1] * wt_cd4_200 +
              get(paste0("prev_200_350_scen",x))[,-1] * wt_cd4_250_350 +
              get(paste0("prev_350_scen",x))[,-1] * wt_cd4_350  ) %>% 
    addYearCol(., horizon_year = horizon_year) 
  
  assign(paste0("hiv_pos_scen",x),hiv_pos)
  
}

##############################################################################

# HIV-NEGATIVE (males and females combined), by CD4+ count and ART status

##############################################################################

#MULTIPLY BY UTILITY WEIGHT

for (x in 1:3) {
  
  hiv_neg <- ((1- (get(paste0("prev_ART_scen",x))[,-1] +    # HIV-negative is 1 - HIV Prevalence
                get(paste0("prev_200_scen",x))[,-1] +
                get(paste0("prev_200_350_scen",x))[,-1] +
                get(paste0("prev_350_scen",x))[,-1] )) 
                      * wt_hiv_neg )  %>%   #multiply by utility weight (in this case = 1)
    addYearCol(., horizon_year = horizon_year)
  
  assign(paste0("hiv_neg_scen",x),hiv_neg)
  
}

# Remove excess DFs

rm(list=ls(pattern="prev"))
rm(hiv_neg, hiv_pos)

##############################################################################

# Annual QALYs - sum HIV positive and HIV negative, then discount

##############################################################################


for (x in 1:3) {
  
  qaly <- ( ( get(paste0("hiv_neg_scen",x))[,-1] + get(paste0("hiv_pos_scen",x))[,-1] )
  
    # MULTIPLY BY POPULATION
    * get(paste0("pop_scen",x))[,-1] ) %>% 
    
    # DISCOUNT
    discount(., discount_rate = discount_rate) %>%  
    
    addYearCol(., horizon_year = horizon_year) %>% 
    recalcFuns(.) # recalculate mean, min, max
  
  assign(paste0("qaly_scen",x),qaly)
  
}

rm(qaly)


# Cumulative 

for (x in 1:3) {
  
  qaly <- get(paste0("qaly_scen",x))[,-1] %>% 
    transmute_at(1:28, ~cumsum(.)) %>% 
    addYearCol(., horizon_year = horizon_year) 
  
  assign(paste0("qaly_cum_scen",x),qaly)
  
}


###################################################################################################################

# Annual raw qalys gained

qaly_gained_raw_scen2 <- (qaly_scen2[,-1] - qaly_scen1[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max 

qaly_gained_raw_scen3 <- (qaly_scen3[,-1] - qaly_scen1[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Annual percent qalys gained

qaly_gained_pct_scen2 <- (((qaly_scen2[,-1] - qaly_scen1[,-1] )/qaly_scen2[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

qaly_gained_pct_scen3 <- (((qaly_scen3[,-1] - qaly_scen1[,-1] )/qaly_scen3[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max 

# Cumulative raw qalys gained

qaly_cum_gained_raw_scen2 <- ( qaly_cum_scen2[,-1] - qaly_cum_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

qaly_cum_gained_raw_scen3 <- ( qaly_cum_scen3[,-1] - qaly_cum_scen1[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative percent qalys gained

qaly_cum_gained_pct_scen2 <- (((qaly_cum_scen2[,-1] - qaly_cum_scen1[,-1] )/qaly_cum_scen2[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

qaly_cum_gained_pct_scen3 <- (((qaly_cum_scen3[,-1] - qaly_cum_scen1[,-1] )/qaly_cum_scen3[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

###################################################################################################################

# Export CSVs

csv_list <- list("qaly_scen1",
                 "qaly_scen2",
                 "qaly_scen3",
                 "qaly_cum_scen1",
                 "qaly_cum_scen2",
                 "qaly_cum_scen3",
                 "qaly_gained_raw_scen2",
                 "qaly_gained_raw_scen3",
                 "qaly_gained_pct_scen2",
                 "qaly_gained_pct_scen3",
                 "qaly_cum_gained_raw_scen2",
                 "qaly_cum_gained_raw_scen3",
                 "qaly_cum_gained_pct_scen2",
                 "qaly_cum_gained_pct_scen3")


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"results/dr3/qaly/",x,".csv"), row.names = F))

# Remove excess DFs, source for next script

rm(list=ls()[! ls() %in% c("main_path","cea_path","helper_path")])
source(paste0(helper_path, "helper.R"))

