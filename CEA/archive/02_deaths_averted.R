
# Compile deaths averted across DO-ART scenarios, combined for Males and Females
# MSahu
# September 10, 2020

# NOTE that annual mortality is calculated per 100K persons
# NOTE: just find and replace from Script 01 : incidence --> mortality, cases --> deaths


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

# IMPORT MORTALITY (males and females combined, age 15-79)

##############################################################################

for (x in 1:3) {
  
  mortality <- read.csv(paste0(main_path,"Scenario",x,"/HIV_mortality_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4)   %>% 
    filter(year>=2020) 
  
  assign(paste0("mortality",x),mortality)
  
}


##############################################################################

# REVISED METHOD: Calculate raw and percentage deaths averted ; discounted

##############################################################################

# New deaths

for (x in 1:3) {
  
  deaths <- ((get(paste0("mortality",x))[,-1]*get(paste0("pop_scen",x))[,-1])/100000)  %>% # Get deaths given mortality rate per 100k  
    
    # DISCOUNT
    discount(., discount_rate = discount_rate) %>%  
    addYearCol(., horizon_year = horizon_year) %>% 
    recalcFuns(.) # recalculate mean, min, max
  
  assign(paste0("deaths_scen",x),deaths)
  
}

# Cumulative deaths

for (x in 1:3) {
  
  deaths_cum <- cumsum(get(paste0("deaths_scen",x))[,-1])  %>% 
    addYearCol(., horizon_year = horizon_year) 
  
  assign(paste0("deaths_cum_scen",x),deaths_cum)
  
}

###################################################################################################################

# Annual raw deaths averted 

deaths_averted_raw_scen2 <- (deaths_scen1[,-1] - deaths_scen2[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

deaths_averted_raw_scen3 <- (deaths_scen1[,-1] - deaths_scen3[,-1] )  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Annual percent deaths averted

deaths_averted_pct_scen2 <- (((deaths_scen1[,-1] - deaths_scen2[,-1] )/deaths_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

deaths_averted_pct_scen3 <- (((deaths_scen1[,-1] - deaths_scen3[,-1] )/deaths_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative raw deaths averted 

deaths_cum_averted_raw_scen2 <- ( deaths_cum_scen1[,-1] - deaths_cum_scen2[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

deaths_cum_averted_raw_scen3 <- ( deaths_cum_scen1[,-1] - deaths_cum_scen3[,-1] ) %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

# Cumulative percent deaths averted

deaths_cum_averted_pct_scen2 <- (((deaths_cum_scen1[,-1] - deaths_cum_scen2[,-1] )/deaths_cum_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

deaths_cum_averted_pct_scen3 <- (((deaths_cum_scen1[,-1] - deaths_cum_scen3[,-1] )/deaths_cum_scen1[,-1] )*100)  %>% 
  addYearCol(., horizon_year = horizon_year) %>% 
  recalcFuns(.) # recalculate mean, min, max

###################################################################################################################

# Export CSVs

csv_list <- list("deaths_scen1",
                 "deaths_scen2",
                 "deaths_scen3",
                 "deaths_cum_scen1",
                 "deaths_cum_scen2",
                 "deaths_cum_scen3",
                 "deaths_averted_raw_scen2",
                 "deaths_averted_raw_scen3",
                 "deaths_averted_pct_scen2",
                 "deaths_averted_pct_scen3",
                 "deaths_cum_averted_raw_scen2",
                 "deaths_cum_averted_raw_scen3",
                 "deaths_cum_averted_pct_scen2",
                 "deaths_cum_averted_pct_scen3")


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"results/",dr,"/deaths/",x,".csv"), row.names = F))

# Remove excess DFs, source for next script

rm(list=ls()[! ls() %in% c("main_path","cea_path","helper_path")])
source(paste0(helper_path, "helper.R"))

