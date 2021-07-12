# ICER table 
# Msahu
# OCt 22, 2020

# ICERs are calculated for each of the 25 parameter sets, and then mean min and max are taken

###########################################################################################

# Incremental Cumulative Costs 

lev <- c("total","pc")

for (i in 2) {
  for (l in lev) {
      
      t <- read.csv(paste0(cea_path,"results/dr3/costs/inc_costs_cum_",l,"_scen",i,".csv"))
      
      x <- t[t$year==horizon_year,-1] %>% # Extract all sets for year of interest, minus year column
      # Recalculate mean, min, max                                                  cv./,
      mutate(mean=rowMeans(.[,4:28]),
             min=apply(.[,4:28], 1,FUN=min), 
             max=apply(.[,4:28], 1,FUN=max)) 
      
      print(paste(horizon_year, "inc costs", l, "scenario", i, ":",round(x$mean,0)))  
      print(paste(horizon_year, "inc costs", l, "scenario MIN", i, ":",round(x$min,0)))  
      print(paste(horizon_year, "inc costs", l, "scenario MAX", i, ":",round(x$max,0)))  
      
      assign(paste0("costsInc_",l,i),x)

  }
}  


# All Effects - UNDISCOUNTED OUTCOMES.

vars <- c("cases", "deaths", "qaly", "daly")
lev <- c("pct", "raw")

for (v in vars) {
  for (i in 2) {
    for (l in lev) {
  
      ifelse(v!="qaly",
        t <- read.csv(paste0(cea_path,"results/dr0/",v,"/",v,"_cum_averted_",l,"_scen",i,".csv")),
        t <- read.csv(paste0(cea_path,"results/dr0/",v,"/",v,"_cum_gained_",l,"_scen",i,".csv")))
      
      x <- t[t$year==horizon_year,-1] %>% # Extract all sets minus year, minus year column
        # Recalculate mean, min, max
        mutate(mean=rowMeans(.[,4:28]),
               min=apply(.[,4:28], 1,FUN=min), 
               max=apply(.[,4:28], 1,FUN=max)) 
      
      print(paste(horizon_year, v, l, "scenario", i, ":",round(x$mean,0)))  
      print(paste(horizon_year, v, l, "scenario MIN", i, ":",round(x$min,0)))  
      print(paste(horizon_year, v, l, "scenario MAX", i, ":",round(x$max,0)))  
      
        assign(paste0(v,"Inc_",l,i),x)

    }
    
    icerDF <- ( get(paste0("costsInc_total",i)) / x ) %>% 
      mutate(year = horizon_year ) %>% 
      select(year, everything())  %>% 
      # Recalculate mean, min, max
      mutate(mean=rowMeans(.[,5:29]),
             min=apply(.[,5:29], 1,FUN=min), 
             max=apply(.[,5:29], 1,FUN=max)) 

    print(paste(horizon_year, v, "ICER scenario", i, ":",round(icerDF$mean,0)))
    print(paste(horizon_year, v, "ICER scenario MIN", i, ":",round(icerDF$min,0)))
    print(paste(horizon_year, v, "ICER scenario MAX", i, ":",round(icerDF$max,0)))
    
  }  
}

rm(t, x)



