#=====================================================================================

# Calculate average annual costs per year
# Author: Mita Sahu
# Last update: July 13, 2021

# NOTE: Must be connected to IHME VPN for currency conversion function 

# To edit any parameters, go to Excel files in "CEA/parameters"

# TO do: 
# 1. SENSITIVITY ANALYSIS
# 2. Total cost - add VMMC?

# =====================================================================================

# Setup ; Read inflation converter from H drive

setwd("/")

library(pacman)

source('H:/repos/fgh/FUNCTIONS/currency_conversion.R')

setwd("C:/Users/msahu/Documents/Other_Research/DO_ART/HHCoM/")


# =======================================================================================

# Inflate literature costs ; assign parameter values

# =======================================================================================

cost_params <- read.csv(paste0(cea_path,"parameters/costs.csv")) %>% 
  mutate(iso3="USA") 

cost_inflated <- currency_conversion(data = cost_params,
                      col.value = "param_val",
                      col.loc = "iso3", 
                      base.year = 2020, 
                      col.currency="currency", 
                      col.currency.year="year",
                      converter.version = 5.2)

costs <- cost_inflated[["param_val"]]
names(costs) <- cost_inflated[["param_name"]]

# Here: add 8 versions of the costs, with lower bounds , etc.

# ========================================================================================

# IMPORT POPULATION

# ========================================================================================

raw_pop_list <- vector(mode = "list")

gender <- c("combined","males","females")

for (v in names(version)) {  
  
  for (g in gender) {
    
    for (x in scenarios) {
      
      temp <- read.csv(paste0(main_path,x,"/PopulationSize_",g,"_aged15-79.csv"), header=F) %>% 
        
        # Set up column names and restrict to 2020:2060
        setNames(df_names) %>% 
        filter(year>=2020)  
      
      raw_pop_list[[length(raw_pop_list) + 1]] <- temp
      
      names(raw_pop_list)[length(raw_pop_list)] <- paste(x, g, v, sep = ".")
    }
  }
}


#===============================================================================================

# IMPORT PREVALENCE - Calculate Raw Prevalent Cases, by CD4+ and ART status

#===============================================================================================

cd4_cat <- c("onART", "CD4_below200_noART", "CD4_200-350_noART", "CD4_350plus_noART")

raw_prev_list <- vector(mode = "list")

for (v in names(version)) {  
  
  for (c in cd4_cat) {
    
    for (x in scenarios) {
      
      if (c != "onART") {  # For people not on ART, don't disaggregate by gender
        
      g <- "combined"
      
      prev <- read.csv(paste0(main_path,x,"/HIV_prevalence_", g, "_aged15-79_",c,".csv"), header=F) %>% 
        
        # Set up column names and restrict to 2020:2060
        setNames(df_names) %>% 
        filter(year>=2020)  %>%  select(-year)
      
      # Total population
      
      pop_size <-  raw_pop_list[[paste(x, g, v, sep =".")]][, -1] 
      
      # Raw prevalent cases = prevalence * population
      
      raw_prev <- prev * pop_size  
      raw_prev <- raw_prev %>% addYearCol() %>% recalcFuns()
      
      raw_prev_list[[length(raw_prev_list) + 1]] <- raw_prev
      
      names(raw_prev_list)[length(raw_prev_list)] <- paste(x, g, c,v, sep = ".")
      
      }
      
      if (c == "onART") {  # For people on ART only, disaggregate by gender (needed for ART costs)
        
        for (g in gender) {
            
            prev <- read.csv(paste0(main_path,x,"/HIV_prevalence_",g,"_aged15-79_",c,".csv"), header=F) %>% 
              
              # Set up column names and restrict to 2020:2060
              setNames(df_names) %>% 
              filter(year>=2020)  %>%  select(-year)
            
            # Total population
            
            pop_size <-  raw_pop_list[[paste(x, g, v, sep =".")]][, -1] 
            
            # Raw prevalent cases = prevalence * population
            
            raw_prev <- prev * pop_size
            raw_prev <- raw_prev %>% addYearCol() %>% recalcFuns()
            
            raw_prev_list[[length(raw_prev_list) + 1]] <- raw_prev
            
            names(raw_prev_list)[length(raw_prev_list)] <- paste(x, g, c,v, sep = ".")
        }
      }
    }
  }
}

rm(prev, pop_size, raw_prev, temp)


#================================================================================================

# HOSPITALIZATION COSTS

# HIV prevalence * average annual cost of hospitalization, by CD4+ count

# Plus incident cases, for 1/2 year - at rate of >350+ untreated

#================================================================================================

temp_list <- vector(mode = "list")
costs_hosp_list <- vector(mode = "list")

inc_cost = costs["hosp_CD4_350plus_noART"]  # Hospitalization cost for incident cases

for (v in names(version)) {  
    
    for (x in scenarios) {
      
      for (d in names(discount_rate)) {
        
        # Prevalent cases (undiscounted)
      
        for (c in cd4_cat) {
      
          # Average annual cost of hospitalization, by CD4+
          
          cost <- costs[paste0("hosp_", c)]
          
          # Prevalent cases, by CD4+
          
          prev <-  raw_prev_list[[paste(x, "combined", c, v, sep =".")]][, -1] 
          
          # Calculate hospitalization costs
          
          hosp_cost = prev * cost 
            
          temp_list[[length(temp_list) + 1]] <- hosp_cost
          
          names(temp_list)[length(temp_list)] <- paste(x, v, c, sep = ".")
          
      }
      
        # Incident cases (undiscounted)
        
        inc <-  outcomes_list[[paste(x, "HIV_incidence", "dr0", v, sep =".")]][, -1] 
        
        hosp_cost = 0.5 * inc * inc_cost 
        
        temp_list[[length(temp_list) + 1]] <- hosp_cost
        
        names(temp_list)[length(temp_list)] <- paste(x, v, "incident_cases", sep = ".")
      
      # Aggregate total costs, by scenario
      
      scenario_hosp_cost <-  Reduce(`+`, temp_list[grep(paste(x,v, sep="."), names(temp_list))] ) 
      
      # Discount 
      
      scenario_hosp_cost <- scenario_hosp_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      costs_hosp_list[[length(costs_hosp_list) + 1]] <- scenario_hosp_cost # Add to list
      
      names(costs_hosp_list)[length(costs_hosp_list)] <- paste(x, "hosp", d, v, sep = ".") # Name
      
      # Reset as empty list
      
      temp_list <- vector(mode = "list")
    }     
  }
}

rm(temp_list, inc_cost, prev, hosp_cost, inc)

#================================================================================================

# TESTING COSTS 

# For Scenario 2 and 2a only, assume one testing campaign every 5 years
# This testing campaign covers HIV negatives plus undiagnosed HIV positives
# Assumes that 75% of people will be reached.

# INPUT: Cara's files on numbers of people tested "Raw_HTC" ; cost and coverage assumptions
# OUTPUT: Total testing cost per year, undiscounted, by scenario

# ??Do we need to add HIV positive not virally suppressed??

#================================================================================================

htc_cost_neg <-  costs["home_testing_pc_neg"]
htc_cost_pos <- costs["home_testing_pc_pos"]

costs_test_list <- vector(mode = "list")

for (v in names(version)) {  
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # POPULATION
      
      if (x %in% c("Scenario1") ) {  # For Scenario #1, zero people on HTC
        
        pop_htc_hivNeg <-  as.data.frame(matrix(0, ncol = 28, nrow = 41)) %>% 
          addYearCol() %>% 
          setNames(df_names) 
        
        pop_htc_hivUndiag <-  as.data.frame(matrix(0, ncol = 28, nrow = 41)) %>% 
          addYearCol() %>% 
          setNames(df_names) 
        
      }
     
      else { # Scenario 2 and 2a
        
        pop_htc_hivNeg <- read.csv(paste0(main_path,x,"/Raw_HTC_combined_hivNeg.csv"), header=F) %>% 
          
          # Set up column names and restrict to 2020:2060
          setNames(df_names) %>% 
          filter(year>=2020)  %>%  select(-year)
        
        pop_htc_hivUndiag <- read.csv(paste0(main_path,x,"/Raw_HTC_combined_hivUndiag.csv"), header=F) %>% 
          
          # Set up column names and restrict to 2020:2060
          setNames(df_names) %>% 
          filter(year>=2020)  %>%  select(-year)
        
      }
      
      # COST
      
      test_cost = pop_htc_hivNeg * htc_cost_neg + pop_htc_hivUndiag * htc_cost_pos
      
      # Discount 
      
      test_cost <- test_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      costs_test_list[[length(costs_test_list) + 1]] <- test_cost # Add to list
      
      names(costs_test_list)[length(costs_test_list)] <- paste(x, "test", d, v, sep = ".") # Name
      
    }     
  }
}

rm(pop_htc_hivUndiag, pop_htc_hivNeg, htc_cost_neg, htc_cost_pos, scenario_hosp_cost, test_cost)


#================================================================================================

# ART COSTS 

# Costs for prevalent cases + costs for cases newly starting ART (for 1/2 year)

# Multiply prevalence of people on ART by the ART costs, using scalar for percentage on community 
# versus clinic-based ART (disaggregated for males and females)

# We use a scalar for percentage of people on ART who are virally suppressed to back-calculate the
# total number of people on ART

# Cases starting / discontinuing ART: calculate costs for 1/2 year (disaggregated by M/F)

#================================================================================================

# COST PARAMATERS  ------------------------------------------------------------------------------

cost_soc_art <- costs["standard_art"]
cost_cb_art <- c(costs["cb_art_y1"], # 1st year
                  rep(costs["cb_art_sub"], 40)) # Subsequent years


# SCALARS  ---------------------------------------------------------------------------------------

# DO ART parameters: assumptions for # people on community versus clinic ART
# See "art_paramaeters" Excel doc for scalars, and "Scenario Details" Word Doc July 8, 2021 for decisions

art_distn <- as.data.frame(readxl::read_excel(paste0(cea_path,"parameters/art_parameters.xls"), sheet = "Table_3")) 

# Also, scalars for percent of peopls on ART who are VS

art_scalars <- as.data.frame(readxl::read_excel(paste0(cea_path,"parameters/art_parameters.xls"), sheet = "Table_2")) 

vs_col <- "Percent who achieve viral suppression of PLHIV who know their status and initiate ART"

FpctVS_soc <- art_scalars[art_scalars$`Treatment intervention`=="Clinic ART" & art_scalars$Gender == "Women", vs_col]

FpctVS_cbArt <- art_scalars[art_scalars$`Treatment intervention`=="Community ART" & art_scalars$Gender == "Women", vs_col]

MpctVS_soc <- art_scalars[art_scalars$`Treatment intervention`=="Clinic ART" & art_scalars$Gender == "Men", vs_col]

MpctVS_cbArt <- art_scalars[art_scalars$`Treatment intervention`=="Community ART" & art_scalars$Gender == "Men", vs_col]

# LOOP ------------------------------------------------------------------------------------------------------
                                   
costs_art_list <- vector(mode = "list")
costs_cbArt_list <- vector(mode = "list")
costs_socArt_list <- vector(mode = "list")

for (v in names(version)) {  
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # POPULATION  ----------------------------------------------------------------------------------------
      # Calculated as: Prevalent cases on ART + raw net initiation (starting ART - discontinuing ART)
      
        # Prevalent cases on ART who are virally suppressed (ART + VS), by gender
        
        Fpop_prev_art <- raw_prev_list[[paste(x, "females", "onART", v, sep =".")]][, -1] 
        
        Mpop_prev_art <- raw_prev_list[[paste(x, "males", "onART", v, sep =".")]][, -1] 
        
        # Raw net ART initiation (ART + VS), by gender  (new inititation minus discontinuation)
        
        Fpop_init_art <- read.csv(paste0(main_path,x,"/Raw_net_ART_init_females_aged15-79.csv"), header=F) %>% 
          setNames(df_names) %>%  # Column names
          filter(year>=2020)  %>%  select(-year) # Restrict to 2020:2060
        
        Mpop_init_art <- read.csv(paste0(main_path,x,"/Raw_net_ART_init_males_aged15-79.csv"), header=F) %>% 
          setNames(df_names) %>%  # Column names
          filter(year>=2020)  %>%  select(-year) # Restrict to 2020:2060
        
        # Set up scalars (Proportion of people on community versus SOC ART)
          
        Fscalar_soc <- art_distn[art_distn$Scenario==x & art_distn$Sex == "Women", "Distribution Clinic ART"]
            
        Fscalar_cbArt <- 1 - Fscalar_soc
          
        Mscalar_soc <- art_distn[art_distn$Scenario==x & art_distn$Sex == "Men", "Distribution Clinic ART"]
          
        Mscalar_cbArt <- 1 - Mscalar_soc
        
        # Total population on ART - backcalculated using scalars for % people on ART virally suppressed
        
        if (vs_scalar == "on") {
        
          pop_soc = (Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_soc  + 
                    (Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art)) / MpctVS_soc 
          
          pop_cb_art = (Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_cbArt + 
                       Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art) / MpctVS_cbArt
        }
        
        if (vs_scalar == "off") {
          
          pop_soc = Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)  + 
                    Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art)
          
          pop_cb_art = Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art) + 
                       Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art)
        }
      
      # COST ------------------------------------------------------------------------------------------------
      
      soc_art_cost = cost_soc_art * pop_soc # multiply pop by soc art scalar
      
      cb_art_cost =  cost_cb_art * pop_cb_art # multiply pop by vector of cb-art costs
      
      art_cost = soc_art_cost + cb_art_cost
      
      # Discount 
      
      art_cost <- art_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      soc_art_cost <- soc_art_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      cb_art_cost <- cb_art_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      costs_art_list[[length(costs_art_list) + 1]] <- art_cost # Add to list
      costs_socArt_list[[length(costs_socArt_list) + 1]] <- soc_art_cost # Add to list
      costs_cbArt_list[[length(costs_cbArt_list) + 1]] <- cb_art_cost # Add to list
      
      names(costs_art_list)[length(costs_art_list)] <- paste(x, "art", d, v, sep = ".") # Name
      names(costs_socArt_list)[length(costs_socArt_list)] <- paste(x, "socArt", d, v, sep = ".") # Name
      names(costs_cbArt_list)[length(costs_cbArt_list)] <- paste(x, "cbArt", d, v, sep = ".") # Name
      
    }     
  }
}

rm(art_scalars, art_distn, vs_col,
   Fscalar_soc, Fscalar_cbArt, Mscalar_soc, Mscalar_cbArt,
   Fpop_init_art, Fpop_prev_art, Mpop_init_art, Mpop_prev_art,
   FpctVS_soc, FpctVS_cbArt, MpctVS_soc, MpctVS_cbArt, 
   pop_soc, pop_cb_art, art_cost, soc_art_cost, cb_art_cost)

#================================================================================================

# VMMC COSTS

# Use files titled "Raw_VMMC_male_ages15-79.csv"

#================================================================================================

costs_vmmc_list <- vector(mode = "list")
vmmc_cost = costs["circumcision"]

for (v in names(version)) {  
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # Population
      
      vmmc_pop <- read.csv(paste0( get(paste0(v, "_path")), x,"/Raw_VMMC_male_ages15-79.csv"), header=F) %>% 
        setNames(df_names) %>%  # Column names
        filter(year>=2020)  %>%  select(-year) # Restrict to 2020:2060
      
      # Cost
      
      cost_vmmc = vmmc_pop * vmmc_cost
      
      # Discount 
      
      cost_vmmc <- cost_vmmc %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      costs_vmmc_list[[length(costs_vmmc_list) + 1]] <- cost_vmmc # Add to list
      
      names(costs_vmmc_list)[length(costs_vmmc_list)] <- paste(x, "vmmc", d, v, sep = ".") # Name
      
    }
  }
}
      
rm(vmmc_pop, cost_vmmc)

#==============================================================================================================

# Calculate Total & Cumulative Costs

# =============================================================================================================

# Total Costs

costs_total_list <- vector(mode = "list")

for (v in names(version)) {  
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
    
      total_cost <-  
        
        costs_hosp_list[[paste(x, "hosp", "dr0",v, sep = ".")]][, -1] +
        
        costs_test_list[[paste(x, "test", "dr0",v, sep = ".")]][, -1] +
        
        costs_art_list[[paste(x, "art", "dr0",v, sep = ".")]][, -1] 
        
        if (vmmc_cost == "on") {
          
          costs_vmmc_list[[paste(x, "vmmc", "dr0",v, sep = ".")]][, -1] 
          
        }
      
      # Discount 
      
      total_cost <- total_cost %>% 
        discount(., discount_rate = discount_rate[d]) %>%
        addYearCol() %>% recalcFuns()
      
      costs_total_list[[length(costs_total_list) + 1]] <- total_cost # Add to list
      
      names(costs_total_list)[length(costs_total_list)] <- paste(x, "total", d, v, sep = ".") # Name
      
    }  
  }
}

# Cumulative Costs

costs_cum_list <- vector(mode = "list")

costs_cum_list <- lapply(costs_total_list, df_cumsum)

#==============================================================================================================

# Incremental Cumulative Costs

# =============================================================================================================

# Set up list

costs_cum_incr_list <- vector(mode = "list")

for (v in names(version)) {  

  for (x in scenarios) {
      
    for (d in names(discount_rate)) {
        
        # Substract : Costs in Alternate Scenario - Costs in Scenario 1
        
        temp <- ( costs_cum_list[[paste(x,"total", d,v, sep =".")]][, -1] -  # Other Scenario 
          
                    costs_cum_list[[paste("Scenario1", "total",d,v, sep =".")]][, -1] ) %>%   # Scenario 1
                    
                  addYearCol() %>%  recalcFuns()
        
        costs_cum_incr_list[[length(costs_cum_incr_list) + 1]] <- temp
        
        names(costs_cum_incr_list)[length(costs_cum_incr_list)] <- paste(x, "total", d, v, sep = ".")
        
    }
  }
}

rm(temp)
