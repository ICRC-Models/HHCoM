
# Compile HIV cases and HIV deaths averted across DO-ART scenarios, combined for Males and Females
# Update: June 24, 2021
# MSahu

# NOTE that annual incidence is calculated per 100 HIV-negative persons 

# To Do:
# Check discounting rates

#=============================================================================================================

# IMPORT Raw Cases & Deaths (males and females combined, age 15-79)

# ============================================================================================================

measure <- c("HIV_incidence", "HIV_mortality")

# Create list of dataframes
# One per scenario (4), discount rate (3), and measure (2)

raw_outcomes_list <- vector(mode = "list")

for (m in measure) {

  for (x in scenarios) {
    
    for (d in discount_rate) {
      
      temp <- read.csv(paste0(main_path, x, "/Raw_", m, "_combined_ages15-79.csv"), header=F) %>% 
        
        # Set up column names and restrict to 2020:2060
        setNames(paste0("s",-3:25)) %>% 
        rename(year=1, mean=2, min=3, max=4) %>% 
        filter(year>=2020) %>% 
        
        # Discount 
        discount(., discount_rate = d) %>% 
        addYearCol()
      
      raw_outcomes_list[[length(raw_outcomes_list) + 1]] <- temp
      
    }
  }
}

names(raw_outcomes_list) <- apply(expand.grid(names(discount_rate), scenarios, measure), 1, paste, collapse=".")

# NEED TO CHECK DISCOUNT RATES worked properly, and everything is named correctly! (use HE book)

#==============================================================================================================

# Calculate Cumulative Raw Cases & Deaths (males and females combined, age 15-79)

# =============================================================================================================

# New list of dataframes - cumulative sum
# One per scenario (4), discount rate (3), and measure (2)

cum_outcomes_list <- vector(mode = "list")

df_cumsum <- function(df) {
 
  df <- cumsum(df[, -1]) %>% 
    addYearCol() %>% 
    recalcFuns()
  
    return(df)
}

cum_outcomes_list <- lapply(raw_outcomes_list, df_cumsum)

#==============================================================================================================

# Calculate Cases & Deaths Averted (males and females combined, age 15-79)

# =============================================================================================================

# Annual -------------------------------------------------------------------------------------------------------

annual_averted_list <- vector(mode = "list")

for (m in measure) {
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # Substract : Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario
      
      temp <- ( raw_outcomes_list[[paste(d,"Scenario1",m, sep =".")]][, -1] -  # Scenario 1
              
              raw_outcomes_list[[paste(d,x,m, sep =".")]][, -1] ) %>%  # Other Scenario
      
              addYearCol() %>%  recalcFuns()
        
      annual_averted_list[[length(annual_averted_list) + 1]] <- temp
      
    }
  }
}

names(annual_averted_list) <- apply(expand.grid(names(discount_rate), scenarios, measure), 1, paste, collapse=".")

# Cumulative -----------------------------------------------------------------------------------------------------

cum_averted_list <- vector(mode = "list")

for (m in measure) {
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # Substract : Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario
      
      temp <- ( cum_outcomes_list[[paste(d,"Scenario1",m, sep =".")]][, -1] -  # Scenario 1
                  
                  cum_outcomes_list[[paste(d,x,m, sep =".")]][, -1] ) %>%  # Other Scenario
        
        addYearCol() %>%  recalcFuns()
      
      cum_averted_list[[length(cum_averted_list) + 1]] <- temp
      
    }
  }
}

names(cum_averted_list) <- apply(expand.grid(names(discount_rate), scenarios, measure), 1, paste, collapse=".")

# Percents (Cumulative) -------------------------------------------------------------------------------------------

pct_cum_averted_list <- vector(mode = "list")

for (m in measure) {
  
  for (x in scenarios) {
    
    for (d in names(discount_rate)) {
      
      # Percent: 100 * (Cases/Deaths in Scenario 1 - Cases/Deaths in Alternate Scenario ) / (Cases/Deaths in Scenario 1)
      
      temp <- ( 100 * (( cum_outcomes_list[[paste(d,"Scenario1",m, sep =".")]][, -1] -  # Scenario 1
                  
                 cum_outcomes_list[[paste(d,x,m, sep =".")]][, -1] ) /  # Other Scenario
                 
                 cum_outcomes_list[[paste(d,"Scenario1",m, sep =".")]][, -1])) %>%  # Scenario 1
        
        addYearCol() %>%  recalcFuns()
      
      pct_cum_averted_list[[length(pct_cum_averted_list) + 1]] <- temp
      
    }
  }
}

names(pct_cum_averted_list) <- apply(expand.grid(names(discount_rate), scenarios, measure), 1, paste, collapse=".")

rm(temp)
