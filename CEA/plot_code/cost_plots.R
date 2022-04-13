##############################################################################################
library(dplyr)
library(ggplot2)
library(ggpubr)

title_size = 15
pres_size = 13.5

lancet_red = "#AD002AFF"
lancet_dblue = "#00468BFF"
lancet_green = "#00A08A"
darjeeling_orange = "#F98400"
darjeeling_yellow = "#F2AD00"
darjeeling_green = "#00A08A"
darjeeling_blue = "#5BBCD6"
darjeeling_red = "#FF0000"
lancet_grey = "#ADB6B6FF"


temp <- c(5,7,6,4,8)
barplot(temp, col="#00A08A")
  
# PLOTS

# Stacked bar chart of breakdown of costs - Scenario 1

cbArt_costs <- costs_cbArt_list$Scenario1.cbArt.dr0.main %>% 
  select(year,mean) %>% 
  rename(`Community ART` = mean)

socArt_costs <- costs_socArt_list$Scenario1.socArt.dr0.main %>% 
  select(year, mean) %>% 
  rename(`Clinic ART` = mean)

hosp_costs <- costs_hosp_list$Scenario1.hosp.dr0.main %>% 
  select(year,mean) %>% 
  rename(Hospitalization=mean)

test_costs <- costs_test_list$Scenario1.test.dr0.main %>% 
  select(year,mean) %>% 
  rename(Testing=mean)

scen1_costs <- socArt_costs %>% 
  left_join(hosp_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("Clinic ART", "Hospitalization"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category`=variable)

scen1_UI <- costs_total_list$Scenario1.total.dr0.main %>% 
  select(-c(mean, min, max)) %>% 
  reshape2::melt(id.vars = "year", measure.vars = c(2:26)) %>% 
  group_by(year) %>% 
  summarize(mean = mean(value),
            min = min(value),
            max = max(value))

p1 <- ggplot() + 
  geom_bar(data = scen1_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year),
           position="stack", stat="identity") +
  geom_errorbar(data = scen1_UI, aes(x = year, ymin = min/1e6, ymax = max/1e6), width = 0.6) +
  scale_fill_manual(values = c(darjeeling_green,
                               lancet_dblue)) +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_blank(),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size),
        legend.position = "bottom")  +
  ylim(0,650)  +
  ggtitle("A) Total Cost for Standard of Care Scenario")


# Stacked bar chart of breakdown of costs - Scenario 2

cbArt_costs <- costs_cbArt_list$Scenario2.cbArt.dr0.main %>% 
  select(year,mean) %>% 
  rename(`Community ART` = mean)

socArt_costs <- costs_socArt_list$Scenario2.socArt.dr0.main %>% 
  select(year, mean) %>% 
  rename(`Clinic ART` = mean)

hosp_costs <- costs_hosp_list$Scenario2.hosp.dr0.main %>% 
  select(year,mean) %>% 
  rename(Hospitalization=mean)

test_costs <- costs_test_list$Scenario2.test.dr0.main %>% 
  select(year,mean) %>% 
  rename(`Home Testing \n Campaign`=mean)

scen2_costs <- cbArt_costs %>% 
  left_join(socArt_costs, by = "year") %>% 
  left_join(hosp_costs, by="year") %>% 
  left_join(test_costs, by="year") %>% 
  reshape2::melt(id.var="year", measure.vars=c("Community ART","Clinic ART", "Hospitalization","Home Testing \n Campaign"), value.name = "Cost") %>% 
  rename(Year=year,
         `Cost Category` = variable)

scen2_UI <- costs_total_list$Scenario2.total.dr0.main %>% 
  select(-c(mean, min, max)) %>% 
  reshape2::melt(id.vars = "year", measure.vars = c(2:26)) %>% 
  group_by(year) %>% 
  summarize(mean = mean(value),
            min = min(value),
            max = max(value))

p3 <- ggplot() + 
  geom_bar(data = scen2_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year),
           position="stack", stat="identity") +
  geom_errorbar(data = scen2_UI, aes(x = year, ymin = min/1e6, ymax = max/1e6), width = 0.6) +
  scale_fill_manual(values = c(darjeeling_yellow,
                               darjeeling_green,
                               lancet_dblue,
                               lancet_red)) +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_blank(),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size),
        legend.position = "bottom")  +
  ylim(0,650) +
  ggtitle("B) Total Cost for Home Testing + Community ART Scenario")

  # Stacked bar chart of breakdown of costs - Scenario 2a
 
  socArt_costs <- costs_socArt_list$Scenario2a.socArt.dr0.main %>% 
    select(year, mean) %>% 
    rename(`Clinic ART` = mean)
  
  hosp_costs <- costs_hosp_list$Scenario2a.hosp.dr0.main %>% 
    select(year,mean) %>% 
    rename(Hospitalization=mean)
  
  test_costs <- costs_test_list$Scenario2a.test.dr0.main %>% 
    select(year,mean) %>% 
    rename(`Home Testing \n Campaign`=mean)
  
  scen2a_costs <- socArt_costs %>% 
    left_join(hosp_costs, by="year") %>% 
    left_join(test_costs, by="year") %>% 
    reshape2::melt(id.var="year", measure.vars=c("Clinic ART", "Hospitalization","Home Testing \n Campaign"), value.name = "Cost") %>% 
    rename(Year=year,
           `Cost Category`=variable)
  
  
  scen2a_UI <- costs_total_list$Scenario2a.total.dr0.main %>% 
    select(-c(mean, min, max)) %>% 
    reshape2::melt(id.vars = "year", measure.vars = c(2:26)) %>% 
    group_by(year) %>% 
    summarize(mean = mean(value),
              min = min(value),
              max = max(value))
  
  p2 <- ggplot() + 
    geom_bar(data = scen2a_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year),
             position="stack", stat="identity") +
    geom_errorbar(data = scen2a_UI, aes(x = year, ymin = min/1e6, ymax = max/1e6), width = 0.6) +
    scale_fill_manual(values = c(darjeeling_green,
                                 lancet_dblue,
                                 lancet_red)) +
    ylab("Total cost (millions, 2020 USD)") + theme_bw() +
    theme(axis.title=element_text(size=pres_size),
          axis.text=element_text(size=pres_size),
          legend.title = element_blank(),
          legend.text = element_text(size = pres_size),
          title = element_text(size = title_size),
          legend.position = "bottom")  +
    ylim(0,650) +
    ggtitle("Total Cost for Home Testing + Clinic-Based Care Scenario")
  

#--------------------------------------------------------------------------------

# Total costs for scenario 1 versus 2

scen1_costs <- costs_total_list$Scenario1.total.dr0.main %>% 
  select(year,mean) %>% 
  rename(`Cost` =mean,
         `Year` = year) %>% 
  mutate(Scenario = "Standard of Care")

scen2_costs <- costs_total_list$Scenario2.total.dr0.main  %>% 
  select(year,mean) %>% 
  rename(`Cost` =mean,
         `Year` = year) %>% 
  mutate(Scenario = "Home Testing + \n Community ART")

total_costs <- scen1_costs %>% 
  rbind(scen2_costs) %>% 
  group_by(Year) %>% mutate(cost_diff = (Cost - lag(Cost, default = first(Cost)))/1e6 )

year_cost_saving = min(total_costs[total_costs$cost_diff<0 & !is.na(total_costs$cost_diff), "Year"])

# Create UIs

scen1_reshape <- costs_total_list$Scenario1.total.dr0.main %>% 
  select(-c(mean, min, max)) %>% 
  reshape2::melt(id.vars = "year", measure.vars = c(2:26)) %>% 
  rename(param_set = variable, cost_scen1 = value)

scen2_reshape <- costs_total_list$Scenario2.total.dr0.main %>% 
  select(-c(mean, min, max)) %>% 
  reshape2::melt(id.vars = "year", measure.vars = c(2:26)) %>% 
  rename(param_set = variable, cost_scen2 = value)

DIFF_UI <- scen1_reshape %>% 
  left_join(scen2_reshape, by = c("year", "param_set")) %>% 
  mutate(diff = cost_scen2 - cost_scen1) %>% 
  group_by(year) %>% 
  summarize(mean_diff = mean(diff),
            min_diff = min(diff),
            max_diff = max(diff))
  
rm(scen1_reshape, scen2_reshape)


p4 <- ggplot(total_costs, aes(y=Cost/1e6, x=Year, fill = Scenario)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(darjeeling_orange, lancet_grey)) +
  annotate(geom="text", x=2049, y=465, label="First year of \n cost savings", color= lancet_red, fontface = "bold") +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  geom_vline(xintercept = year_cost_saving, linetype = "dashed", color = lancet_red, size = 1.2) +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size)) +
  ylim(0, 500) +
  ggtitle("Standard of Care Scenario versus \n Home Testing + Community ART Scenario")

p4 <- ggplot(total_costs, aes(y=Cost/1e6, x=Year, fill = Scenario)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(darjeeling_orange, lancet_grey)) +
  annotate(geom="text", x=2049, y=465, label="First year of \n cost savings", color= lancet_red, fontface = "bold") +
  ylab("Total cost (millions, 2020 USD)") + theme_bw() +
  geom_vline(xintercept = year_cost_saving, linetype = "dashed", color = lancet_red, size = 1.2) +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size)) +
  ylim(0, 500) +
  ggtitle("Standard of Care Scenario versus \n Home Testing + Community ART Scenario")

p4_revised <- ggplot() + 
  geom_bar(data = total_costs, aes(y=cost_diff, x=Year), position="dodge", stat="identity", color = lancet_grey, fill = darjeeling_orange) +
  geom_errorbar(data = DIFF_UI, aes(x = year, ymin = min_diff/1e6, ymax = max_diff/1e6), width = 0.6) + 
  annotate(geom="text", x=2055.5, y=67, label="First year of projected \n cost savings = 2049", color= lancet_red, fontface = "bold", size = 4) +
  ylab("Incremental cost (millions, 2020 USD)") + theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = year_cost_saving + 0.5, linetype = "dashed", color = lancet_red, size = 1.2) +
  theme(axis.title=element_text(size=pres_size),
        axis.text=element_text(size=pres_size),
        legend.title = element_text(size = pres_size),
        legend.text = element_text(size = pres_size),
        title = element_text(size = title_size)) +
 # scale_y_continuous(limits = c(-5, 80), breaks = seq(-10, 80, 10)) +
  ggtitle("C) Cost Difference for Home Testing + Community ART \n compared with Standard of Care Scenario")


# fig2 <- ggarrange(p1, p2, p3, p4,  ncol = 1,
#                   label.x = 1, label.y = 1)
# 
# ggsave(plot = fig2, file = paste0(dir, "CEA/figures/Figure2.png"), device = "png",
#        width = 8, height = 12)
# 
# 
# fig2_v2 <- ggarrange(p1, p3,p4,  ncol = 1,
#                   label.x = 1, label.y = 1)
# 
# ggsave(plot = fig2_v2, file = paste0(dir, "CEA/figures/Figure2_v2.png"), device = "png",
#        width = 8, height = 12)

fig2_v3 <- ggarrange(p1, p3,p4_revised,  ncol = 1,
                     label.x = 1, label.y = 1)

ggsave(plot = fig2_v3, file = paste0(dir, "CEA/figures/Figure2.pdf"), device = "pdf",
       width = 8, height = 12)

ggsave(plot = fig2_v3, file = paste0(dir, "CEA/figures/Figure2.eps"), device = "eps",
       width = 8, height = 12)


#--------------------------------------------------------------------------------

# INCLUDING VMMC COSTS

#--------------------------------------------------------------------------------
# 
# 
# # Stacked bar chart of breakdown of costs - Scenario 1
# 
# cbArt_costs <- costs_cbArt_list$Scenario1.cbArt.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(`Community ART` = mean)
# 
# socArt_costs <- costs_socArt_list$Scenario1.socArt.dr0.vmmc %>%
#   select(year, mean) %>%
#   rename(`Clinic ART` = mean)
# 
# hosp_costs <- costs_hosp_list$Scenario1.hosp.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(Hospital=mean)
# 
# test_costs <- costs_test_list$Scenario1.test.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(Testing=mean)
# 
# vmmc_costs <- costs_vmmc_list$Scenario1.vmmc.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(VMMC=mean)
# 
# scen1_costs <- cbArt_costs %>%
#   left_join(socArt_costs, by = "year") %>%
#   left_join(hosp_costs, by="year") %>%
#   left_join(test_costs, by="year") %>%
#   left_join(vmmc_costs, by="year") %>%
#   reshape2::melt(id.var="year", measure.vars=c("Community ART","Clinic ART", "Hospital","Testing", "VMMC"), value.name = "Cost") %>%
#   rename(Year=year,
#          `Cost Category`=variable)
# 
# ggplot(scen1_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year)) +
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(values = c("darkorange1","#42B540FF", "#00468BFF", "#AD002AFF", "cyan4")) +
#   ylab("Total cost (millions, 2020 USD)") + theme_bw() +
#   theme(axis.title=element_text(size=pres_size),
#         axis.text=element_text(size=pres_size),
#         legend.title = element_text(size = pres_size),
#         legend.text = element_text(size = pres_size),
#         title = element_text(size = title_size))  +
#   ylim(0,500)  +
#   ggtitle("Standard of Care \n Scenario")
# 
# # Stacked bar chart of breakdown of costs - Scenario 2
# 
# cbArt_costs <- costs_cbArt_list$Scenario2.cbArt.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(`Community ART` = mean)
# 
# socArt_costs <- costs_socArt_list$Scenario2.socArt.dr0.vmmc %>%
#   select(year, mean) %>%
#   rename(`Clinic ART` = mean)
# 
# hosp_costs <- costs_hosp_list$Scenario2.hosp.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(Hospital=mean)
# 
# test_costs <- costs_test_list$Scenario2.test.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(Testing=mean)
# 
# vmmc_costs <- costs_vmmc_list$Scenario2.vmmc.dr0.vmmc %>%
#   select(year,mean) %>%
#   rename(VMMC=mean)
# 
# scen2_costs <- cbArt_costs %>%
#   left_join(socArt_costs, by = "year") %>%
#   left_join(hosp_costs, by="year") %>%
#   left_join(test_costs, by="year") %>%
#   left_join(vmmc_costs, by="year") %>%
#   reshape2::melt(id.var="year", measure.vars=c("Community ART","Clinic ART", "Hospital","Testing", "VMMC"), value.name = "Cost") %>%
#   rename(Year=year,
#          `Cost Category`=variable)
# 
# ggplot(scen2_costs, aes(fill=`Cost Category`, y=Cost/1e6, x=Year)) +
#   geom_bar(position="stack", stat="identity") +
#  scale_fill_manual(values = c("darkorange1","#42B540FF", "#00468BFF", "#AD002AFF", "cyan4")) +
#   ylab("Total cost (millions, 2020 USD)") + theme_bw() +
#   theme(axis.title=element_text(size=pres_size),
#         axis.text=element_text(size=pres_size),
#         legend.title = element_text(size = pres_size),
#         legend.text = element_text(size = pres_size),
#         title = element_text(size = title_size))  +
#   ylim(0,500) +
#   ggtitle("HTC + Community ART \n Scenario")
# 
