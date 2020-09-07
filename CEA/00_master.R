# MASTER SCRIPT
# September 7, 2020
#MSahu

# Set up 

rm(list=ls())

library("dplyr")

main_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/DoArtOutputs/"
cea_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/CEA/"

# Set up helper discount functions

discount_rate <- .03   
discount <- function(year_discount,discount_rate) {
  1/((1 + discount_rate)^year_discount)
}
discounter <- function(x,discount_amt) {
  x*discount_amt
}

# NEXT STEPS:
# Males/females
# Vary Discount Rate 0-5%