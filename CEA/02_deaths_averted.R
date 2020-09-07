# Code to compile Deaths Averted across DoART scenarios, combined Males and Females
# MSahu
# September 7, 2020

# NOTE that annual mortality is calculated per 100K persons
# NOTE: just find and replace from Script 01 : mortality --> mortality, cases --> deahts

source("00_master.R")
library("dplyr")

# Import mortality - combined across gender

for (x in 1:3) {
  
  mortality <- read.csv(paste0(main_path,"Scenario",x,"/HIV_mortality_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) 
  
  assign(paste0("mortality",x),mortality)
  
}

# Check mean, min, max are correct for all datasets -- ALL OK

mortality_check <- mortality1 %>% mutate(mean2=rowMeans(.[,5:29]),
                                         min2=apply(.[,5:29],1,FUN=min), 
                                         max2=apply(.[,5:29],1,FUN=max))
print(all.equal(mortality_check$mean, mortality_check$mean2))
print(all.equal(mortality_check$min, mortality_check$min2))
print(all.equal(mortality_check$max, mortality_check$max2))

mortality_check <- mortality2 %>% mutate(mean2=rowMeans(.[,5:29]),
                                         min2=apply(.[,5:29],1,FUN=min),   
                                         max2=apply(.[,5:29],1,FUN=max))
print(all.equal(mortality_check$mean,mortality_check$mean2))
print(all.equal(mortality_check$min,mortality_check$min2))
print(all.equal(mortality_check$max,mortality_check$max2))

mortality_check <- mortality3 %>% mutate(mean2=rowMeans(.[,5:29]),
                                         min2=apply(.[,5:29],1,FUN=min),   
                                         max2=apply(.[,5:29],1,FUN=max))
print(all.equal(mortality_check$mean,mortality_check$mean2))
print(all.equal(mortality_check$min,mortality_check$min2))
print(all.equal(mortality_check$max,mortality_check$max2))

rm(mortality_check)


# DISCOUNT mortality
# We discount effects in order to be consistent with discounting of costs (p. 110 in Drummond et all)
# We discount 3% here. Need to vary in sensitivity analysis from 0% to 5%

DFnames <- c("mortality1","mortality2","mortality3")

for (x in DFnames) {
  
  mortality_discounted <- get(x) %>% 
    filter(year>=2020) %>% 
    mutate(year_discount=0:40,   # set 2020 to Year 0
           discount_amt=discount(discount_rate,year_discount)) %>% 
    mutate_at(vars(c(mean,min,max,5:29)),~discounter(.,discount_amt))
  
  assign(paste0(x,"_disc"),mortality_discounted)
  
}

# Calculate deaths averted

deaths_averted_scen2 <- (mortality1_disc[,2:29]-mortality2_disc[,2:29]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) 

write.csv(deaths_averted_scen2, paste0(cea_path,"effects/deaths_averted_scen2.csv"))

deaths_averted_scen3 <- (mortality1_disc[,2:29]-mortality3_disc[,2:29]) %>% 
  mutate(year=2020:2060) %>% 
  select(year,everything()) 

write.csv(deaths_averted_scen3, paste0(cea_path,"effects/deaths_averted_scen3.csv"))


# Remove excess DFs

rm(mortality,mortality1,mortality2,mortality3,mortality_discounted, mortality1_disc,mortality2_disc,mortality3_disc)

