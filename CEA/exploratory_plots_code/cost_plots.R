##############################################################################################

pres_size <- 28

# PLOTS

cea_path <- "C:/Users/msahu/Documents/Other_Research/DO_ART/Code/HHCoM/CEA/"

# Stacked bar chart of breakdown of costs - Scenario 1

art_costs <- read.csv(paste0(cea_path,"costs/ART_scen1.csv")) %>% 
  select(year,mean) %>% 
  rename(ART=mean)

hosp_costs <- read.csv(paste0(cea_path,"costs/hosp_scen1.csv")) %>% 
  select(year,mean) %>% 
  rename(Hospital=mean)

test_costs <- read.csv(paste0(cea_path,"costs/test_scen1.csv")) %>% 
  select(year,mean) %>% 
  rename(Testing=mean)

scen1_costs <- art_costs %>% 
  left_join(hosp_costs, by="year") %>% 
  left_join(test_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("ART","Hospital","Testing"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category`=variable)

ggplot(scen1_costs, aes(fill=`Cost Category`, y=Cost, x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#42B540FF","#00468BFF", "#AD002AFF")) +
  ylab("Cost per capita (2020 USD)") +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size)) +
  ylim(0,60)



# Stacked bar chart of breakdown of costs - Scenario 2

art_costs <- read.csv(paste0(cea_path,"costs/ART_scen2.csv")) %>% 
  select(year,mean) %>% 
  rename(ART=mean) 

hosp_costs <- read.csv(paste0(cea_path,"costs/hosp_scen2.csv")) %>% 
  select(year,mean) %>% 
  rename(Hospital=mean)

test_costs <- read.csv(paste0(cea_path,"costs/test_scen2.csv")) %>% 
  select(year,mean) %>% 
  rename(Testing=mean)

scen2_costs <- art_costs %>% 
  left_join(hosp_costs, by="year") %>% 
  left_join(test_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("ART","Hospital","Testing"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category`=variable)

ggplot(scen2_costs, aes(fill=`Cost Category`, y=Cost, x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#42B540FF","#00468BFF", "#AD002AFF")) +
  ylab("Cost per capita (2020 USD)") +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size)) 



# Annual costs (per capita, full population)

costs <- read.csv(paste0(cea_path,"costs/cost_annual_pc_scen1.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual costs per capita (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,100)

costs <- read.csv(paste0(cea_path,"costs/cost_annual_pc_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual costs per capita (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,100)


# Annual costs (total)

costs <- read.csv(paste0(cea_path,"costs/cost_annual_total_scen1.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs/1000000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual costs per capita (in millions of 2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,500)

costs <- read.csv(paste0(cea_path,"costs/cost_annual_total_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs/1000000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual costs per capita (in millions of 2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,500)


# Incremental costs (per capita)

costs <- read.csv(paste0(cea_path,"costs/inc_costs_annual_pc_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual incremental costs per capita (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  geom_hline(yintercept = 0,linetype="dashed",color="darkgrey")

costs <- read.csv(paste0(cea_path,"costs/inc_costs_cum_pc_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Incremental cumulative costs per capita (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))  +
  geom_hline(yintercept = 0,linetype="dashed",color="darkgrey")

# Incremental costs (total)

costs <- read.csv(paste0(cea_path,"costs/inc_costs_annual_total_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs/1000000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual incremental costs (millions of 2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  geom_hline(yintercept = 0,linetype="dashed",color="darkgrey") +
  scale_y_continuous(labels = function(x) format(x, scientific = F)) 

costs <- read.csv(paste0(cea_path,"costs/inc_costs_cum_pc_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative costs per capita (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))  +
  geom_hline(yintercept = 0,linetype="dashed",color="darkgrey")


costs <- read.csv(paste0(cea_path,"costs/inc_costs_cum_total_scen2.csv")) %>% 
  select(-X) %>% 
  # reshape long
  melt(id="year",value.name="costs") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=costs, aes(x=year, y=costs,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative costs - total population (2020 USD)") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18))  +
  geom_hline(yintercept = 0,linetype="dashed",color="darkgrey")