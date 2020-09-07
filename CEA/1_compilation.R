# Code to compile DoART scenarios
# MSahu
# July 27, 2020

# Load packages

library("dplyr")
library("readxl")

# Setup

model_outputs <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/DoArtOutputs/"

# Generate total mortality over time per capita

mortality_1 <- read_excel(paste0(model_outputs,"Scenario1/HIV_incidence_females_aged15-79.xlsx"), 
                          sheet="HIVincF", col_names=F) %>% 
  mutate(year=[[1]])
         , 2=mean, 3=min, 4=max)