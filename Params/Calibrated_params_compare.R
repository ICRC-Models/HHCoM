library(haven)
library(readxl)
library(tidyverse)
library(sjlabelled)
library(lubridate)
library(WriteXLS)

init <- read.table("Z:/Kenya_model_optim_calib/HHCoM/Params/paramSets_patternSrch_16DEC19_0.dat")
bounds <- read_excel("Z:/Kenya_model_optim_calib/HHCoM/Params/paramSets_patternSrch_16Dec19_bounds.xlsx")
calib <- read.table("Z:/Kenya_model_optim_calib/HHCoM/Params/negSumLogL_patternSrch_30DEC19_0_params.dat")


View(compare_calib)

init <- init %>% rename(intial = V1) %>% add_column(ID = 1:8)
bounds <- bounds %>% add_column(ID=1:8)
compare_calib <- calib %>% rename(calibrated = V1) %>% add_column(ID = 1:8) %>%
     left_join(init, by = "ID") %>%
     left_join(bounds, by = "ID") %>%
     select(-ID)

View(compare_calib)
