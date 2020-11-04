# ICER and  CEA planes ; time plots
# September 10, 2020
# MSahu

library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(reshape2)

############################################################################################

# CALCULATE COST PER CASE AVERTED, DEATH AVERTED, QALY GAINED

mergeCols <- c("year", "set_name")

# Scenario 2 - incremental

inc_costs2 <- read.csv(paste0(cea_path,"costs/inc_costs_cum_total_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="cost_diff") %>%
  rename(set_name=variable) %>%
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

cases_cum_averted_2 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen2.csv")) %>%
 select(-X) %>%
  reshape2::melt(id="year",value.name="cases_averted") %>%
rename(set_name=variable)

deaths_cum_averted_2 <- read.csv(paste0(cea_path,"effects/deaths/deaths_cum_averted_raw_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="deaths_averted") %>%
  rename(set_name=variable)

qalys_cum_gained_2 <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_gained_raw_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="qaly_gained") %>%
  rename(set_name=variable)

dalys_cum_averted_2 <- read.csv(paste0(cea_path,"effects/dalys/daly_cum_averted_raw_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="daly_averted") %>%
  rename(set_name=variable)

cost_case2_60 <- inc_costs2 %>% 
  left_join(cases_cum_averted_2, by=mergeCols) %>%
  left_join(deaths_cum_averted_2, by=mergeCols) %>% 
  left_join(qalys_cum_gained_2, by=mergeCols) %>% 
  left_join(dalys_cum_averted_2, by=mergeCols) %>% 
   filter(year==2060) %>%
   mutate(scenario="2") %>% 
   mutate(cost_per_case = cost_diff/cases_averted,
         cost_per_death = cost_diff/deaths_averted,
         cost_per_qaly = cost_diff/qaly_gained,
         cost_per_daly = cost_diff/daly_averted)


# Scenario 3 - incremental

inc_costs2 <- read.csv(paste0(cea_path,"costs/inc_costs_cum_total_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="cost_diff") %>%
  rename(set_name=variable) %>%
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

cases_cum_averted_2 <- read.csv(paste0(cea_path,"effects/cases/cases_cum_averted_raw_scen2.csv")) %>%
  select(-X) %>%
  reshape2::melt(id="year",value.name="cases_averted") %>%
  rename(set_name=variable)

cost_case2_60 <- inc_costs2 %>% left_join(cases_cum_averted_2) %>%
  filter(year==2060) %>%
  mutate(scenario="2") %>% 
  mutate(cost_per_case = cost_diff/cases_averted)




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

dalys_cum_1 <- read.csv(paste0(cea_path,"effects/dalys/daly_cum_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) 

cea1 <- costs_cum_1 %>% 
  left_join(cases_cum_1, by=mergeCols) %>% 
  left_join(deaths_cum_1, by=mergeCols) %>% 
  left_join(qalys_cum_1, by=mergeCols) %>% 
  left_join(dalys_cum_1, by=mergeCols) %>% 
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

dalys_cum_2 <- read.csv(paste0(cea_path,"effects/dalys/daly_cum_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) 

cea2 <- costs_cum_2 %>% 
  left_join(cases_cum_2, by=mergeCols) %>% 
  left_join(deaths_cum_2, by=mergeCols) %>% 
  left_join(qalys_cum_2, by=mergeCols) %>% 
  left_join(dalys_cum_2, by=mergeCols) %>% 
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

dalys_cum_3 <- read.csv(paste0(cea_path,"effects/dalys/daly_cum_scen3.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) 

cea3 <- costs_cum_3 %>% 
  left_join(cases_cum_3, by=mergeCols) %>% 
  left_join(deaths_cum_3, by=mergeCols) %>% 
  left_join(qalys_cum_3, by=mergeCols) %>% 
  left_join(dalys_cum_3, by=mergeCols) %>% 
  mutate(scenario="3")

# Merge

rm(cases_cum_1, cases_cum_2, cases_cum_3,
   deaths_cum_1, deaths_cum_2, deaths_cum_3,
   qalys_cum_1, qalys_cum_2, qalys_cum_3,
   dalys_cum_1, dalys_cum_2, dalys_cum_3,
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

# Graph with lines

cost60 <- cea60 %>% 
  reshape2::dcast(set_name~scenario,value.var="cost") %>% 
  rename(cost1="1",
         cost2="2")

case60 <- cea60 %>% 
  reshape2::dcast(set_name~scenario,value.var="cases") %>% 
  rename(cases1="1",
         cases2="2")

death60 <- cea60 %>% 
  reshape2::dcast(set_name~scenario,value.var="deaths") %>% 
  rename(deaths1="1",
         deaths2="2")

qaly60 <- cea60 %>% 
  reshape2::dcast(set_name~scenario,value.var="qalys") %>% 
  rename(qaly1="1",
         qaly2="2")

daly60 <- cea60 %>% 
  reshape2::dcast(set_name~scenario,value.var="dalys") %>% 
  rename(daly1="1",
         daly2="2")

cea60_lines <- cost60 %>% 
  left_join(case60, by="set_name") %>% 
  left_join(death60, by="set_name") %>% 
  left_join(qaly60, by="set_name") %>% 
  left_join(daly60, by="set_name") %>% 
  mutate(size=ifelse(set_name=="mean","mean", "other set"))

rm(cost60, case60, death60, qaly60, daly60)


# 2060 Plots with Lines 

ggplot(data=cea60_lines) +
  geom_point(aes(x=cases1/1000, y=cost1/1000000, size=size), color="blue") +
  geom_point(aes(x=cases2/1000, y=cost2/1000000, size=size), color="red") +
  geom_segment(aes(x = cases1/1000, xend = cases2/1000, y = cost1/1000000, yend = cost2/1000000), color="black") +
  scale_size_manual(values=c(5,2)) +
  xlab("Cases (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea60_lines) +
  geom_point(aes(x=deaths1/1000, y=cost1/1000000, size=size), color="blue") +
  geom_point(aes(x=deaths2/1000, y=cost2/1000000, size=size), color="red") +
  geom_segment(aes(x = deaths1/1000, xend = deaths2/1000, y = cost1/1000000, yend = cost2/1000000), color="black") +
  scale_size_manual(values=c(5,2)) +
  xlab("Deaths (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 


ggplot(data=cea60_lines) +
  geom_point(aes(x=qaly1/1000, y=cost1/1000000, size=size), color="blue") +
  geom_point(aes(x=qaly2/1000, y=cost2/1000000, size=size), color="red") +
  geom_segment(aes(x = qaly1/1000, xend = qaly2/1000, y = cost1/1000000, yend = cost2/1000000), color="black") +
  scale_size_manual(values=c(5,2)) +
  xlab("QALYs (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggplot(data=cea60_lines) +
  geom_point(aes(x=daly1/1000, y=cost1/1000000, size=size), color="blue") +
  geom_point(aes(x=daly2/1000, y=cost2/1000000, size=size), color="red") +
  geom_segment(aes(x = daly1/1000, xend = daly2/1000, y = cost1/1000000, yend = cost2/1000000), color="black") +
  scale_size_manual(values=c(5,2)) +
  xlab("DALYs (in thousands) ") +
  ylab("Costs (in millions of USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 


####################################################################################################

# Net Health Benefits (QALYS)

threshold <- 750

costs_cum_inc_2 <- read.csv(paste0(cea_path,"costs/inc_costs_cum_total_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="costs") %>% 
  rename(set_name=variable)

qalys_cum_gained_2 <- read.csv(paste0(cea_path,"effects/qaly/qaly_cum_gained_raw_scen2.csv")) %>% 
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable)

dalys_cum_averted_2 <- read.csv(paste0(cea_path,"effects/daly/daly_cum_averted_raw_scen2.csv")) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable)


NHB_qaly <- costs_cum_inc_2  %>%  left_join(qalys_cum_gained_2, by=c("year","set_name")) %>% 
  mutate(nhb = qalys - (costs/threshold)) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=NHB_qaly, aes(x=year, y=nhb, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net Health Benefits - QALYs") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

NHB_daly <- costs_cum_inc_2  %>%  left_join(dalys_cum_averted_2, by=c("year","set_name")) %>% 
  mutate(nhb = (dalys* (-1)) + (costs/threshold)) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

NHB_daly_REVISED <- costs_cum_inc_2  %>%  left_join(dalys_cum_averted_2, by=c("year","set_name")) %>% 
  mutate(nhb = (dalys - (costs/threshold))) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=NHB_daly_REVISED, aes(x=year, y=nhb, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net Health Benefits - DALYs") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))

# Net Monetary Benefits (DALYS)


threshold <- 750


NMB_qaly <- costs_cum_inc_2  %>%  left_join(qalys_cum_gained_2, by=c("year","set_name")) %>% 
  mutate(nmb = (qalys*threshold) + costs ) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))


ggplot(data=NMB_qaly, aes(x=year, y=nmb, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net Monetary Benefits - QALYs") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 



NMB_daly <- costs_cum_inc_2  %>%  left_join(dalys_cum_averted_2, by=c("year","set_name")) %>% 
  mutate(nmb = ((dalys*(-1))*threshold) + costs ) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))


NMB_daly_REVISED <- costs_cum_inc_2  %>%  left_join(dalys_cum_averted_2, by=c("year","set_name")) %>% 
  mutate(nmb = ((dalys*threshold) - costs )) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))


ggplot(data=NMB_daly_REVISED, aes(x=year, y=nmb, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net Monetary Benefits - DALYs") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 


# Net DALY burden per scenario


threshold <- 750

costs_cum_1 <- read.csv(paste0(cea_path,"costs/cost_cum_total_scen1.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="costs") %>% 
  rename(set_name=variable)

daly_cum_1 <- read.csv(paste0(cea_path,"effects/daly/daly_cum_scen1.csv")) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable)


net_daly_1 <- costs_cum_1 %>%  left_join(daly_cum_1, by=c("year","set_name")) %>% 
  mutate(daly_burden = dalys + (costs/threshold) ) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=net_daly_1, aes(x=year, y=daly_burden, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net DALY Burden, Scenario 1") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 



costs_cum_2 <- read.csv(paste0(cea_path,"costs/cost_cum_total_scen2.csv")) %>% 
  select(-X) %>% 
  reshape2::melt(id="year",value.name="costs") %>% 
  rename(set_name=variable)

daly_cum_2 <- read.csv(paste0(cea_path,"effects/daly/daly_cum_scen2.csv")) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable)


net_daly_2 <- costs_cum_2 %>%  left_join(daly_cum_2, by=c("year","set_name")) %>% 
  mutate(daly_burden = dalys + (costs/threshold) ) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=net_daly_2, aes(x=year, y=daly_burden, group=set_name)) +
  geom_line(aes(color=set_name)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Net DALY Burden, Scenario 2") +
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
