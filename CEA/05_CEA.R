# ICER and  CEA planes ; time plots
# September 10, 2020
# MSahu

library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(reshape2)

# # Scenario 2 - incremental
# 
# inc_costs2 <- read.csv(paste0(cea_path,"costs/inc_costs_annual_total_scen2.csv")) %>% 
#   select(-X) %>% 
#   reshape2::melt(id="year",value.name="cost_diff") %>% 
#   rename(set_name=variable) %>% 
#   mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))
# 
# cases_cum_averted_2 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen2.csv")) %>% 
#   select(-X) %>% 
#   reshape2::melt(id="year",value.name="cases_averted") %>% 
#   rename(set_name=variable) 
# 
# cost_case2_60 <- inc_costs2 %>% left_join(cases_cum_averted_2) %>% 
#   filter(year==2060) %>% 
#   mutate(scenario="2")
# 
# # Scenario 3 - incremental 
# 
# 
# inc_costs3 <- read.csv(paste0(cea_path,"costs/inc_costs_annual_total_scen3.csv")) %>% 
#   select(-X) %>% 
#   reshape2::melt(id="year",value.name="cost_diff") %>% 
#   rename(set_name=variable) %>% 
#   mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))
# 
# cases_cum_averted_3 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen3.csv")) %>% 
#   select(-X) %>% 
#   reshape2::melt(id="year",value.name="cases_averted") %>% 
#   rename(set_name=variable) 
# 
# cost_case3_60 <- inc_costs3 %>% left_join(cases_cum_averted_3) %>% 
#   filter(year==2060) %>% 
#   mutate(scenario="3")
# 
# 
# # Merge
# 
# cost_case60 <- cost_case2_60 %>%  full_join(cost_case3_60)
# 
# 
# ggplot(data=cost_case60, aes(x=cases_averted, y=cost_diff,group=scenario)) +
#   geom_point(aes(color=scenario,size=size)) +
#   scale_size_manual(values=c(5,2)) +
#   xlab("Incremental effectiveness (cases averted)") +
#   ylab("Incremental costs") +
#   theme_cowplot() +
#   theme(axis.title=element_text(size=18),
#         axis.text=element_text(size=18)) 
# 
# 
# cost_case30 <- costs %>% left_join(cases_averted_raw)  %>% filter(year==2030)


########################################################################################

# COMBINED DATASET  - CUMULATIVE

mergeCols <- c("year","set_name")

# Scenario 1 

costs_cum_1 <- read.csv(paste0(cea_path,"costs/cost_cum_total_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cost") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name=="mean","mean value","other sets")) 

cases_cum_1 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cases") %>% 
  rename(set_name=variable) 

deaths_cum_1 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="deaths") %>% 
  rename(set_name=variable) 

qalys_cum_1 <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) 

cea1 <- costs_cum_1 %>% 
  left_join(cases_cum_1, by=mergeCols) %>% 
  left_join(deaths_cum_1, by=mergeCols) %>% 
  left_join(qalys_cum_1, by=mergeCols) %>% 
  mutate(scenario="1") 

# Scenario 2

costs_cum_2 <- read.csv(paste0(cea_path,"costs/cost_cum_total_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cost") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name=="mean","mean value","other sets")) 

cases_cum_2 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cases") %>% 
  rename(set_name=variable) 

deaths_cum_2 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="deaths") %>% 
  rename(set_name=variable) 

qalys_cum_2 <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) 

cea2 <- costs_cum_2 %>% 
  left_join(cases_cum_2, by=mergeCols) %>% 
  left_join(deaths_cum_2, by=mergeCols) %>% 
  left_join(qalys_cum_2, by=mergeCols) %>% 
  mutate(scenario="2")

# Scenario 3 - PLACEHOLDER ONLY

costs_cum_3 <- read.csv(paste0(cea_path,"costs/cost_cum_total_scen3.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cost") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name=="mean","mean value","other sets")) 

cases_cum_3 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_scen3.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="cases") %>% 
  rename(set_name=variable) 

deaths_cum_3 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_scen3.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="deaths") %>% 
  rename(set_name=variable) 

qalys_cum_3 <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_scen3.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) 

cea3 <- costs_cum_3 %>% 
  left_join(cases_cum_3, by=mergeCols) %>% 
  left_join(deaths_cum_3, by=mergeCols) %>% 
  left_join(qalys_cum_3, by=mergeCols) %>% 
  mutate(scenario="3")

# Merge

rm(cases_cum_1, cases_cum_2, cases_cum_3,
   deaths_cum_1, deaths_cum_2, deaths_cum_3,
   qalys_cum_1, qalys_cum_2, qalys_cum_3,
   costs_cum_1, costs_cum_2, costs_cum_3)

#2060 plots

cea60 <- cea1 %>% 
  full_join(cea2) %>% 
  filter(year==2060)

ggplot(data=cea60, aes(x=cases/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Cases (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea60, aes(x=deaths/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Deaths (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea60, aes(x=qalys/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("QALYs (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

#2030 plots

cea30 <- cea1 %>% 
  full_join(cea2) %>% 
  filter(year==2030)

ggplot(data=cea30, aes(x=cases/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Cases (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea30, aes(x=deaths/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Deaths (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea30, aes(x=qalys/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("QALYs (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 


#2045 plots

cea45 <- cea1 %>% 
  full_join(cea2) %>% 
  filter(year==2030)

ggplot(data=cea45, aes(x=cases/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Cases (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea45, aes(x=deaths/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("Deaths (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea45, aes(x=qalys/1000, y=cost/1000000,group=scenario)) +
  geom_point(aes(color=scenario,size=size)) +
  scale_size_manual(values=c(5,2)) +
  xlab("QALYs (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

####################################################################################################

# Annual Percent Deaths Averted 

deaths <- read.csv(paste0(cea_path,"effects/deaths/deaths_averted_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_deaths_averted_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=deaths, aes(x=year, y=pct_deaths_averted_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage deaths averted (annual)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,50)

deaths_cum <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_averted_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_deaths_averted_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=deaths_cum, aes(x=year, y=pct_deaths_averted_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage deaths averted (cumulative)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))  +
  ylim(0,50)





# Annual QALYs gained 

qalys <- read.csv(paste0(cea_path,"effects/qalys/qaly_gained_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_qalys_gained_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys, aes(x=year, y=pct_qalys_gained_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage QALYs gained (annual)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

qalys_cum <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_gained_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_qalys_gained_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys_cum, aes(x=year, y=pct_qalys_gained_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage QALYs gained (cumulative)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 
