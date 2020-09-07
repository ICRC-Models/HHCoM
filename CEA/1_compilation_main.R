# Code to compile DoART scenarios
# MSahu
# July 27, 2020

# Load packages

rm(list=ls())

library("dplyr")
library("readxl")

# Setup

model_outputs <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/DoArtOutputs/"

# Generate total cumulative mortality over time per capita for each Scenario

# COMBINED (later disaggregate by sex)

for (x in 1:3) {

m <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_mortality_combined_aged15-79.xlsx"), 
                          sheet="HIVmortC", col_names=F) %>% 
  rename(year=1, mean=2, min=3, max=4) %>% 
  mutate(cum_mean=cumsum(mean), 
         cum_min=cumsum(min),
         cum_sum=cumsum(max)) %>% 
  select(year,grep("cum",colnames(.)))

assign(paste0("m",x),m)

}

# Merge
mortality <- m1 %>%  
  left_join(m2, by="year") %>%  
  left_join(m3, by="year") %>% 
  #hand clean
  rename(cum_mean1 = cum_mean.x,
         cum_min1 = cum_min.x,
         cum_max1 = cum_sum.x,
         cum_mean2 = cum_mean.y,
         cum_min2 = cum_min.y,
         cum_max2 = cum_sum.y,
         cum_mean3 = cum_mean,
         cum_min3 = cum_min,
         cum_max3 = cum_sum) %>% 
  mutate(diff_2_1_mean = cum_mean2-cum_mean1,
         diff_2_1_min = cum_min2-cum_min1,
         diff_2_1_max = cum_max2-cum_max1,
         diff_3_1_mean = cum_mean3-cum_mean1,
         diff_3_1_min = cum_min3-cum_min1,
         diff_3_1_max = cum_max3-cum_max1) %>% 
  select(year,grep("diff",colnames(.)))

rm(m,m1,m2,m3)


# Generate total cumulative incidence over time per capita for each Scenario

for (x in 1:3) {
  
  i <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_incidence_combined_aged15-79.xlsx"), 
                  sheet="HIVincC", col_names=F) %>% 
    rename(year=1, mean=2, min=3, max=4) %>% 
    mutate(cum_mean=cumsum(mean), 
           cum_min=cumsum(min),
           cum_sum=cumsum(max)) %>% 
    select(year,grep("cum",colnames(.)))
  
  assign(paste0("i",x),i)
  
}

# Merge
incidence <- i1 %>%  
  left_join(i2, by="year") %>%  
  left_join(i3, by="year") %>% 
  #hand clean
  rename(cum_mean1 = cum_mean.x,
         cum_min1 = cum_min.x,
         cum_max1 = cum_sum.x,
         cum_mean2 = cum_mean.y,
         cum_min2 = cum_min.y,
         cum_max2 = cum_sum.y,
         cum_mean3 = cum_mean,
         cum_min3 = cum_min,
         cum_max3 = cum_sum) %>% 
  mutate(diff_2_1_mean = cum_mean2-cum_mean1,
         diff_2_1_min = cum_min2-cum_min1,
         diff_2_1_max = cum_max2-cum_max1,
         diff_3_1_mean = cum_mean3-cum_mean1,
         diff_3_1_min = cum_min3-cum_min1,
         diff_3_1_max = cum_max3-cum_max1) %>% 
  select(year,grep("diff",colnames(.)))


rm(i,i1,i2,i3)



# Prevalence

# on ART


for (x in 1:3) {
  
  m <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_prevalence_combined_aged15-79_onART.xlsx"), 
                  sheet="HIVprevC", col_names=F) %>% 
    rename(year=1, mean=2, min=3, max=4)  
  
  assign(paste0("m",x),m)
  
}

prev_ART <- m1 %>%  
  left_join(m2, by="year") %>%  
  left_join(m3, by="year") %>% 
  #hand clean
  rename(mean1 = mean.x,
         min1 = min.x,
         max1 = max.x,
         mean2 = mean.y,
         min2 = min.y,
         max2 = max.y,
         mean3 = mean,
         min3 = min,
         max3 = max) 

rm(m,m1,m2,m3)


m <- read.csv(paste0(model_outputs,"Scenario1/CD4_dist_initART_males_aged15-79CD4 200-350.csv"))

#CD4 < 200

for (x in 1:3) {
  
  m <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4 below200 noART.xlsx"), 
                  sheet="HIVprevC", col_names=F) %>% 
    rename(year=1, mean=2, min=3, max=4)  
  
  assign(paste0("m",x),m)
  
}

prev_200 <- m1 %>%  
  left_join(m2, by="year") %>%  
  left_join(m3, by="year") %>% 
  #hand clean
  rename(mean1 = mean.x,
         min1 = min.x,
         max1 = max.x,
         mean2 = mean.y,
         min2 = min.y,
         max2 = max.y,
         mean3 = mean,
         min3 = min,
         max3 = max) 

rm(m,m1,m2,m3)


# CD4 200-350


for (x in 1:3) {
  
  m <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4 200-350 noART.xlsx"), 
                  sheet="HIVprevC", col_names=F) %>% 
    rename(year=1, mean=2, min=3, max=4)  
  
  assign(paste0("m",x),m)
  
}

prev_200_350 <- m1 %>%  
  left_join(m2, by="year") %>%  
  left_join(m3, by="year") %>% 
  #hand clean
  rename(mean1 = mean.x,
         min1 = min.x,
         max1 = max.x,
         mean2 = mean.y,
         min2 = min.y,
         max2 = max.y,
         mean3 = mean,
         min3 = min,
         max3 = max) 

rm(m,m1,m2,m3)

# CD4 350+


for (x in 1:3) {
  
  m <- read_excel(paste0(model_outputs,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4 350plus noART.xlsx"), 
                  sheet="HIVprevC", col_names=F) %>% 
    rename(year=1, mean=2, min=3, max=4)  
  
  assign(paste0("m",x),m)
  
}

prev_350 <- m1 %>%  
  left_join(m2, by="year") %>%  
  left_join(m3, by="year") %>% 
  #hand clean
  rename(mean1 = mean.x,
         min1 = min.x,
         max1 = max.x,
         mean2 = mean.y,
         min2 = min.y,
         max2 = max.y,
         mean3 = mean,
         min3 = min,
         max3 = max) 

rm(m,m1,m2,m3)
