# ICER and  CEA planes ; time plots
# September 10, 2020
# MSahu

library(ggplot2)
library(reshape2)
library(cowplot)


# Annual Percent Cases Averted 

cases <- read.csv(paste0(cea_path,"effects/cases/cases_averted_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_cases_averted_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=cases, aes(x=year, y=pct_cases_averted_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage cases averted (annual)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,50)

cases_cum <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_pct_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="pct_cases_averted_2") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=cases_cum, aes(x=year, y=pct_cases_averted_2,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Percentage cases averted (cumulative)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))  +
  ylim(0,50)


# Cost per case averted

cases_averted_raw <- read.csv(paste0(cea_path,"effects/cases/cases_averted_raw_scen2.csv")) %>% 
  select(-X) %>% 
  melt(id="year",value.name="cases_averted_2") %>% 
  rename(set_name=variable) 

cost_case <- costs %>% left_join(cases_averted_raw) %>% 
  mutate(cost_case=cost_diff/cases_averted_2)

ggplot(data=cost_case, aes(x=year, y=cost_case,group=set_name)) +
  geom_point(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Incremental Average Annual Cost / Annual Cases Averted") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 


cost_case30 <- costs %>% left_join(cases_averted_raw)  %>% filter(year==2030)





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
