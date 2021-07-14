##############################################################################################
library(dplyr)
library(ggplot2)
library(ggpubr)

title_size = 50
pres_size = 40

# PLOTS

# Stacked bar chart of breakdown of costs - Scenario 1

cbArt_costs <- costs_cbArt_list$Scenario1.cbArt.dr3.main %>% 
  select(year,mean) %>% 
  rename(`Community ART` = mean)

socArt_costs <- costs_socArt_list$Scenario1.socArt.dr3.main %>% 
  select(year, mean) %>% 
  rename(`Clinic ART` = mean)

hosp_costs <- costs_hosp_list$Scenario1.hosp.dr3.main %>% 
  select(year,mean) %>% 
  rename(Hospital=mean)

test_costs <- costs_test_list$Scenario1.test.dr3.main %>% 
  select(year,mean) %>% 
  rename(Testing=mean)

scen1_costs <- cbArt_costs %>% 
  left_join(socArt_costs, by = "year") %>% 
  left_join(hosp_costs, by="year") %>% 
  left_join(test_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("Community ART","Clinic ART", "Hospital","Testing"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category`=variable)

ggplot(scen1_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("darkorange1","#42B540FF", "#00468BFF", "#AD002AFF")) +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size))  +
  ylim(0,500)  +
  ggtitle("Standard of Care \n Scenario")

# Stacked bar chart of breakdown of costs - Scenario 2

cbArt_costs <- costs_cbArt_list$Scenario2.cbArt.dr3.main %>% 
  select(year,mean) %>% 
  rename(`Community ART` = mean)

socArt_costs <- costs_socArt_list$Scenario2.socArt.dr3.main %>% 
  select(year, mean) %>% 
  rename(`Clinic ART` = mean)

hosp_costs <- costs_hosp_list$Scenario2.hosp.dr3.main %>% 
  select(year,mean) %>% 
  rename(Hospital=mean)

test_costs <- costs_test_list$Scenario2.test.dr3.main %>% 
  select(year,mean) %>% 
  rename(Testing=mean)

scen2_costs <- cbArt_costs %>% 
  left_join(socArt_costs, by = "year") %>% 
  left_join(hosp_costs, by="year") %>% 
  left_join(test_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("Community ART","Clinic ART", "Hospital", "Testing"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category`=variable)

ggplot(scen2_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("darkorange1","#42B540FF", "#00468BFF", "#AD002AFF")) +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size))  +
  ylim(0,500)  +
  ggtitle("HTC + Community ART \n Scenario")


#--------------------------------------------------------------------------------

# Total costs for scenario 1 versus 2

scen1_costs <- costs_total_list$Scenario1.total.dr3.main %>% 
  select(year,mean) %>% 
  rename(`Cost` =mean,
         `Year` = year) %>% 
  mutate(Scenario = "Standard of Care")

scen2_costs <- costs_total_list$Scenario2.total.dr3.main  %>% 
  select(year,mean) %>% 
  rename(`Cost` =mean,
         `Year` = year) %>% 
  mutate(Scenario = "Home Testing + \n Community ART")

total_costs <- scen1_costs %>% 
  rbind(scen2_costs)

ggplot(total_costs, aes(y=Cost/1e6, x=Year, fill = Scenario)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("blue4","dodgerblue3")) +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size)) +
  ylim(0, 500)
  