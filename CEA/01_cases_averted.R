# Code to compile Cases Averted across DoART scenarios, combined Males and Females
# MSahu
# September 7, 2020

# NOTE that annual incidence is calculated per 100 HIV-negative persons 

# Setup

source("00_master.R")
library("dplyr")

# Import incidence - combined across gender

for (x in 1:3) {

incidence <- read.csv(paste0(main_path,"Scenario",x,"/HIV_incidence_combined_aged15-79.csv"), header=F) %>% 
  setNames(paste0("s",-3:25)) %>% 
  rename(year=1,
         mean=2,
         min=3,
         max=4) 

assign(paste0("incidence",x),incidence)

}

# Check mean, min, max are correct for all datasets -- ALL OK

incidence_check <- incidence1 %>% mutate(mean2=rowMeans(.[,5:29]),
                          min2=apply(.[,5:29],1,FUN=min),   # 1 means it's being applied across rows
                          max2=apply(.[,5:29],1,FUN=max))
print(all.equal(incidence_check$mean,incidence_check$mean2))
print(all.equal(incidence_check$min,incidence_check$min2))
print(all.equal(incidence_check$max,incidence_check$max2))

incidence_check <- incidence2 %>% mutate(mean2=rowMeans(.[,5:29]),
                                          min2=apply(.[,5:29],1,FUN=min),   
                                          max2=apply(.[,5:29],1,FUN=max))
print(all.equal(incidence_check$mean,incidence_check$mean2))
print(all.equal(incidence_check$min,incidence_check$min2))
print(all.equal(incidence_check$max,incidence_check$max2))

incidence_check <- incidence3 %>% mutate(mean2=rowMeans(.[,5:29]),
                                          min2=apply(.[,5:29],1,FUN=min),   
                                          max2=apply(.[,5:29],1,FUN=max))
print(all.equal(incidence_check$mean,incidence_check$mean2))
print(all.equal(incidence_check$min,incidence_check$min2))
print(all.equal(incidence_check$max,incidence_check$max2))

rm(incidence_check)


# DISCOUNT INCIDENCE
# We discount effects in order to be consistent with discounting of costs (p. 110 in Drummond et all)

DFnames <- c("incidence1","incidence2","incidence3")

for (x in DFnames) {

incidence_discounted <- get(x) %>% 
    filter(year>=2020) %>% 
    mutate(year_discount=0:40,   # set 2020 to Year 0
           discount_amt=discount(discount_rate,year_discount)) %>% 
    mutate_at(vars(c(mean,min,max,5:29)),~discounter(.,discount_amt))

assign(paste0(x,"_disc"),incidence_discounted)

}

# Calculate cases averted

cases_averted_scen2 <- (incidence1_disc[,2:29]-incidence2_disc[,2:29]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) 

write.csv(cases_averted_scen2, paste0(cea_path,"effects/cases_averted_scen2.csv"))

cases_averted_scen3 <- (incidence1_disc[,2:29]-incidence3_disc[,2:29]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) 

write.csv(cases_averted_scen3, paste0(cea_path,"effects/cases_averted_scen3.csv"))

# Remove excess DFs

rm(incidence,incidence1,incidence2,incidence3,incidence_discounted, incidence1_disc,incidence2_disc,incidence3_disc)
