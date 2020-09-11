# Checks 
# MSahu
# Sep 10, 2020

# CHECK COLUMNS NAMED AS EXPECTED


###################################################################################################################


for (x in 1:3) {
  
  incidence <- read.csv(paste0(main_path,"Scenario",x,"/HIV_incidence_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4)   %>% 
    filter(year>=2020) 
  
  assign(paste0("incidence",x),incidence)
  
}

# Check mean, min, max are Columns 1-3 as expected for all datasets -- ALL OK

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

########################################################################################################

for (x in 1:3) {
  
  mortality <- read.csv(paste0(main_path,"Scenario",x,"/HIV_mortality_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    rename(year=1,
           mean=2,
           min=3,
           max=4) 
  
  assign(paste0("mortality",x),mortality)
  
}

# Check mean, min, max are Columns 1-3 as expected for all datasets -- ALL OK

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
