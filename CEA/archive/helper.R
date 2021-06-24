# Helper functions and parameters
# November 18, 2020

#############################################################################

## PARAMETERS ##

# Set time horizon

horizon_year <- 2060

# Discount rate of 3%

discount_rate <- .0
dr <- "dr0"


##############################################################################

# Set up DO ART parameters: enrollment assumptions

females_ART_scen1 <- .6223
females_ART_scen2 <- .657
females_ART_scen3 <- .857

males_ART_scen1 <- .3978
males_ART_scen2 <- .655
males_ART_scen3 <- .857

males_enrolment_clinic <- 51/72
males_enrolment_cbART <- 21/72

females_enrolment_clinic <- 70/73
females_enrolment_cbART <- 3/73

# Set up DO ART % tested scalar

DOARTpct_tested <- 0.9

##############################################################################

# FUNCTIONS #

# Add year column from 2020 to time horizon

addYearCol <- function(x, horizon_year) {
  x %>% 
    mutate(year= 2020:horizon_year) %>% 
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


