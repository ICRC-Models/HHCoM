
# Compile HIV cases and HIV deaths averted across DO-ART scenarios, combined for Males and Females
# Update: June 24, 2021
# MSahu

#=============================================================================================================

# IMPORT Raw Cases & Deaths (males and females combined, age 15-79)

# Remember:
# Annual incidence is calculated per 100 HIV-negative persons 
# Annual mortality is calculated per capita
# But use raw cases here

# ============================================================================================================

measure <- c("HIV_incidence", "HIV_mortality")

# Create list of dataframes
# One per scenario (4), discount rate (3), and measure (2)

outcomes_list <- vector(mode = "list")

for (v in names(version)) {  

  for (m in measure) {  
  
    for (x in scenarios) {
      
      for (d in names(discount_rate)) {
        
        temp <- read.csv(paste0(main_path, x, "/Raw_", m, "_combined_ages15-79.csv"), header=F) %>% 
          
          # Set up column names and restrict to 2020:2060
          setNames(df_names) %>% 
          filter(year>=2020) %>% 
          
          # Discount 
          discount(., discount_rate = discount_rate[d]) %>% 
          addYearCol()
        
        outcomes_list[[length(outcomes_list) + 1]] <- temp
        
        names(outcomes_list)[length(outcomes_list)] <- paste(x, m, d, v, sep = ".")
        
      }
    }
  }
}

#==============================================================================================================

# Calculate Cumulative Raw Cases & Deaths (males and females combined, age 15-79)

# =============================================================================================================

# New list of dataframes - cumulative sum
# One per scenario (4), discount rate (3), and measure (2)

outcomes_cum_list <- vector(mode = "list")

outcomes_cum_list <- lapply(outcomes_list, df_cumsum)

#==============================================================================================================

# Calculate Cases & Deaths Averted (males and females combined, age 15-79), and Percents

# =============================================================================================================

# Set up lists

outcomes_avrt_ann_list <- vector(mode = "list")
outcomes_avrt_cum_list <- vector(mode = "list")

for (v in names(version)) {  

  for (m in measure) {
    
    for (x in scenarios) {
      
      for (d in names(discount_rate)) {
        
        # Annual --------------------------------------------------------------------------------------------------
        
        # Substract : Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario
        
        temp <- ( outcomes_list[[paste("Scenario1",m,d,v, sep =".")]][, -1] -  # Scenario 1
                
                outcomes_list[[paste(x,m,d,v, sep =".")]][, -1] ) %>%  # Other Scenario
        
                addYearCol() %>%  recalcFuns()
          
        outcomes_avrt_ann_list[[length(outcomes_avrt_ann_list) + 1]] <- temp
        
        # Cumuluative ---------------------------------------------------------------------------------------------
        
        # Substract : Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario
        
        temp <- ( outcomes_cum_list[[paste("Scenario1",m,d,v, sep =".")]][, -1] -  # Scenario 1
                     
                     outcomes_cum_list[[paste(x,m,d,v, sep =".")]][, -1] ) %>%  # Other Scenario
          
          addYearCol() %>%  recalcFuns()
        
        outcomes_avrt_cum_list[[length(outcomes_avrt_cum_list) + 1]] <- temp
        
        names(outcomes_avrt_cum_list)[length(outcomes_avrt_cum_list)] <- paste(x, m, d, v, sep = ".")
        
      }
    }
  }
}


# Percents (Cumulative) -------------------------------------------------------------------------------------------

outcomes_avrt_cum_pct_list <- vector(mode = "list")

for (v in names(version)) {  

  for (m in measure) {
    
    for (x in scenarios) {
      
      for (d in names(discount_rate)) {
        
        # Percent: 100 * (Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario ) / (Cases/Deaths in Scenario 1)
        
        temp <- ( 100 * (( outcomes_cum_list[[paste("Scenario1",m,d,v, sep =".")]][, -1] -  # Scenario 1
                    
                   outcomes_cum_list[[paste(x,m,d,v, sep =".")]][, -1] ) /  # Other Scenario
                   
                   outcomes_cum_list[[paste("Scenario1",m,d,v, sep =".")]][, -1])) %>%  # Scenario 1
          
          addYearCol() %>%  recalcFuns()
        
        outcomes_avrt_cum_pct_list[[length(outcomes_avrt_cum_pct_list) + 1]] <- temp
        
        names(outcomes_avrt_cum_pct_list)[length(outcomes_avrt_cum_pct_list)] <- paste(x, m, d, v, sep = ".")
        
      }
    }
  }
}

rm(temp)


