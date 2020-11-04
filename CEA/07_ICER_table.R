# ICER table 
# Msahu
# OCt 22, 2020

# TO DO: Fix ICER to calculate across all sets

###########################################################################################

horizon_year <- 2060

# Incremental Cumulative Costs

lev <- c("total","pc")

for (i in 2:3) {
  for (l in lev) {
      
      t <- read.csv(paste0(cea_path,"costs/inc_costs_cum_",l,"_scen",i,".csv"))
      
      x <- t[t$year==horizon_year, "mean"] # Extract mean value for year of interest
      
      assign(paste0("costsInc_",l,i),x)
      print(paste(horizon_year, "costs", l, "scenario", i, ":",round(x,1)))  # Print
  }
}  


# All Effects

vars <- c("cases", "deaths", "qaly", "daly")
lev <- c("pct", "raw")

for (v in vars) {
  for (i in 2:3) {
    for (l in lev) {
  
      ifelse(v!="qaly",
        t <- read.csv(paste0(cea_path,"effects/",v,"/",v,"_cum_averted_",l,"_scen",i,".csv")),
        t <- read.csv(paste0(cea_path,"effects/",v,"/",v,"_cum_gained_",l,"_scen",i,".csv")))
  
      x <- t[t$year==horizon_year, "mean"] # Extract mean value for year of interest
      assign(paste0(v,"Inc_",l,i),x)
      print(paste(horizon_year, v, l, "scenario", i, ":",round(x,1)))  
    }
    
    # NEED TO FIX ("ICERS" - should be calculating the mean across all sets)
    
    icer <- get(paste0("costsInc_total",i)) / x
    print(paste(horizon_year, v, "ICER scenario", i, ":",round(icer,1)))
    assign(paste0(v,"ICER",i),icer)
  }  
}

rm(t, x, icer)





