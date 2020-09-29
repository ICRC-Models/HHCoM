# ICER

#MSahu
#July 28, 2020

library(dplyr)
library(ggplot2)
library(reshape2)

costs2 <- read.csv(paste0(cea_path,"costs/inc_costs_annual_total_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  mutate(scenario=2)


cases2 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="cases_averted") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=2)


deaths2 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_averted_raw_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="deaths_averted") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=2)


qalys2 <- read.csv(paste0(cea_path,"effects/qalys/qaly_gained_raw_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="qalys_gained") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=2)
 


costs3 <- read.csv(paste0(cea_path,"costs/inc_costs_annual_total_scen3.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=3)



cases3 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen3.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="cases_averted") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=3)



deaths3 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_averted_raw_scen3.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="deaths_averted") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=3)


qalys3 <- read.csv(paste0(cea_path,"effects/qalys/qaly_gained_raw_scen3.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="qalys_gained") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal")) %>% 
  filter(year==2060) %>% 
  select(-size) %>% 
  mutate(scenario=3)



mergeCols <- c("year","set_name","scenario")
merged <- costs2 %>% full_join(cases2, by=mergeCols) %>% 
  full_join(deaths2, by=mergeCols) %>% 
  full_join(qalys2, by=mergeCols) %>% 
  full_join(costs3, by=mergeCols) %>% 
  full_join(cases3, by=mergeCols) %>% 
  full_join(deaths3, by=mergeCols) %>% 
  full_join(qalys3, by=mergeCols) %>% 
  mutate(costs= ifelse(scenario==2, costs2,costs3))  %>% 
  mutate(cases_averted= ifelse(scenario==2, cases_averted.x,cases_averted.y)) %>% 
  mutate(deaths_averted= ifelse(scenario==2, deaths_averted.x,deaths_averted.y))                      


icer <- mortality %>% rename(mortality_2_1_diff = diff_2_1_mean) %>% select(year, mortality_2_1_diff) %>% 
  left_join(incidence, by ="year") %>%  rename(incidence_2_1_diff = diff_2_1_mean) %>% 
  select(year, incidence_2_1_diff, mortality_2_1_diff) %>% 
  left_join(cost_diff) %>% 
  filter(year>=2016) %>% 
  mutate(cost_death_averted= cost_diff_1_2/(mortality_2_1_diff*100000)) %>% 
  mutate(cost_case_averted=cost_diff_1_2/(incidence_2_1_diff*100))
  
ggplot(data=merged, aes(x=costs.x, y=cases_averted.x)) +
  geom_point() +
  xlab("Incremental costs (2020 USD)") +
  ylab("Cases Averted") +
  theme_bw()

ggplot(data=icer, aes(x=year, y=cost_case_averted)) +
  geom_line() +
  xlab("Year") +
  ylab("Cost per case averted") +
  theme_bw()


