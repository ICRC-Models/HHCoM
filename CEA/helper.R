# Helper functions 
# Last update: July 1, 2021
# Msahu

# ===================================================================================

# FUNCTIONS #

# Add year column from 2020 to 2060

addYearCol <- function(x) {
  x %>% 
    mutate(year= 2020:2060) %>% 
    select(year, everything()) 
}

# Recalculate the mean, min, and max across 25 parameter sets named s1...s25

recalcFuns <- function(x) {
  x %>% 
    mutate(mean = rowMeans(.[ , grepl( "s\\d+" , names(.) ) ]),
           min = apply(.[ , grepl( "s\\d+" , names(.) ) ] , MARGIN = 1, FUN = min),
           max = apply(.[ , grepl( "s\\d+" , names(.) ) ] , MARGIN = 1, FUN = max))
}

# Discounter

discount <- function(x, discount_rate) {
  x %>% 
    mutate( year_discount = 0:(nrow(.)-1),
            discount_pct = (1/((1 + discount_rate)^year_discount) )) %>% 
    transmute_if( grepl( pattern = "s\\d+", x = names(.) ) , 
                  ~.*discount_pct) %>% 
    recalcFuns(.) %>% 
    select(mean, min, max,everything())
}

# Cumulative sum

df_cumsum <- function(df) {
  
  df <- cumsum(df[, -1]) %>% 
    addYearCol() %>% 
    recalcFuns()
  
  return(df)
}