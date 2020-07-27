reactive_surveillance = function(epi_input_cease, thresh_end, thresh_start, yrs_active_off, yrs_react){
  
  # NOTE: Around September 2019, the epi team at Warwick took care of 
  # simulating the reactive surveillance within the transmission model. 
  # I leave this function for convenience.
  
  ### Inputs
  # epi_input_cease: epi_input array (from the epi team) limited to these outcomes 
  # "Active_S1"    "Active_S2"    "Passive_S1"   "Passive_S2".
  # the dimensions are: iterations, years (years of intervention only), outcome, intervention
  
  # thresh_end: number of cases that would warrant cutting active surveillance
  # thresh_start: number of cases that must be found to trigger AS again (default is 1), 
  # indicate Inf to indicate that reactive surveillance is off
  
  # yrs_active_off: number of years with thresh_end cases to shut down active surveillance
  # yrs_react: number of years reactive surveillance runs before shutting down again
  
  # Dependencies: fcn_end_yr (a function we created)

  ### Note for future: some active surveillance finds
  # might trigger active surveillance in village of origin... 
  # Reactive surveillance might not be exhaustive (smaller population) - how to code that?!?
  # Reactive surveillance might have startup costs or different 
  # costs to active surveillance - for now keep it simple
  # int_units_all - if enhanced passive surveillance is aborted, then some of 
  # the units of the passive surveillance components will change too.
  
  ### Values:
  # 1. as_status: each year: as, off, ras1, ras2 ... # then you can ask for 
  # for the following summaries (perhaps a post-prediction function for it)
  
  # end date #1 (very last year if it never ends)
  # reactive surveillance year #1
  # end date to RAS #1
  # reactive surveillance year #2 
  # end date to RAS #2
  # number of reactive surveillances that had to be done
  
  # 2. epi_input_cease1, counts_int_cease, cost_int_as_cease: each year
  # keep costs by type, by component of intervention, mgmt vs unit too.
  
  ### Intermediate arrays to make:
  # 1. tmp_int_units: from counts_int number of people where each kind of surveillance is administered 
  ## dimensions: iterations, component, interventions, years, type_act 
  ## (components related to AS both main and secondary activities to be modified; 
  ## potentially enhanced passive as well in the future)
  # 2. tmp_cas_all: from epi_input_cease array with cases found through both kinds of surveillance (active and passive)
  ## dimensions: iterations, years, intervention
  # 3. tmp_cas_pass: from epi_input_cease array with cases found through passive surveillance only
  ## dimensions: iterations, years, intervention 
  
  ### Steps in the function: 
  # a. make tmp objects
  # b. could use the rle function to find the year when active surveillance first shuts down
  # c. fill in the first part of as_status with "as" using that year.
  # d. then, starting on that year (when AS shuts down) examine if it should 
  # stay shut down or reactivated (record in as_status)
  # e. what happens to those missed? LATER; for now they all die precisely the year when they were missed.
  
  # a. make tmp objects
  tmp_cas_all = apply(epi_input_cease, c(1,2,4), "sum", na.rm=T)
  tmp_cas_pass = apply(epi_input_cease[,,c("Passive_S1", "Passive_S2"),], c(1,2,4), "sum", na.rm=T)
  
  # b. could use the rle function to find the year when active surveillance first shuts down
  tmp_yrs_abort = sapply(1:dim(tmp_cas_all)[3], 
                         function(d){apply(tmp_cas_all[,,d], 1, "end_yr", yrs_active_off, thresh_end)})
  
  # c. Fill in as_status with as either "on" up until the first year of shutdown.
  as_status = array('as', dim=dim(epi_input_cease)[c(1,2,4)],
                    dimnames = dimnames(epi_input_cease)[c(1,2,4)])
  # for (it in 1:dim(epi_input_cease)[1]){
  #   for (int in 1:dim(epi_input_cease)[3]){
  #     as_status[it,  1:(which(dimnames(epi_input_cease)$years==as.character(tmp_yrs_abort[it,int]))-1), int] = "as"
  #   }
  # }
  
  # d. then, starting on that year (when AS shuts down) examine if it should 
  # stay shut down or reactivated (record in as_status)
  # b=Sys.time()
  for (it in 1:dim(epi_input_cease)[1]){
    for (int in 1:dim(epi_input_cease)[3]){
      
      yr = tmp_yrs_abort[it,int]
      if(yr<max(as.numeric(dimnames(epi_input_cease)$years))){
        as_mode="off"
        # AS situation #1: the off years
        while(yr<=max(as.numeric(dimnames(epi_input_cease)$years))){ 
            while(yr<=max(as.numeric(dimnames(epi_input_cease)$years)) && 
              tmp_cas_pass[it,as.character(yr),int]<thresh_start && as_mode=="off"){
              as_status[it, as.character(yr), int] = "off"
              yr=yr+1
            }
      
        # AS situation 2a: the beginning of the ras years
        # for the first yrs_react years (years of reactive surveillance 
        # before considering shutting down again), it's reactive surveillance
        
      if(yr<=max(as.numeric(dimnames(epi_input_cease)$years))){ # if yr is still under max
          if(thresh_start <= tmp_cas_pass[it,as.character(yr),int]){
            st_yr = min(yr+yrs_react, max(as.numeric(dimnames(epi_input_cease)$years)))
            
            as_status[it, as.character(yr), int] = "off" # the year that the case is found in PS, keep AS off
            
            if (yr+1<=st_yr){
            as_status[it, as.character((yr+1):st_yr), int] = "ras"
            }
            yr = st_yr
            as_mode="ras"
          }
        # after that, you look at the previous three years of tmp_cas_all NOT tmp_cas_pass
        
        # AS situation 2b: the ras years continue (if cases don't go back down)
          if(any(tmp_cas_all[it,as.character((yr-yrs_react):(yr-1)),int]>0) & as_mode=="ras"){
            as_status[it, as.character(yr), int] = "ras"
            yr=yr+1
          } else if (!any(tmp_cas_all[it,as.character((yr-yrs_react):(yr-1)),int]>0) & as_mode=="ras") { 
            # then it goes back to off, so write down the status as off, walk one more time in the yr vector, and turn as_mode off
            as_status[it, as.character(yr), int] = "off"
            yr=yr+1
            as_mode="off"
            # then it will go back to the first while loop
          }
        }
      }
      }  
      
      }
  }
  # a=Sys.time() 11 seconds
  
  return(as_status)
  
  # e. what happens to those missed? LATER.
  # uncaught cases: S2's that eventually make it to passive surveillance, 
  # S1s that go to be S2s and make it to passive surveillance
  # As of 4 June 2019: they all die. Add them up.
  
} # end of reactive_surveillance function





