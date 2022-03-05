# ICER and  CEA planes ; time plots
# MSahu
# Update: 3/4/2022

library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(reshape2)

#========================================================================================================

#  Import data

costs_cum_inc_2 <- costs_cum_incr_list$Scenario2.total.dr3.main %>% 
  select(-c(mean, min, max)) %>% 
  reshape2::melt(id.var = "year") %>% 
  rename(inc_cost = value) %>% 
  mutate(param_set = as.numeric(substring(variable, 2))) %>% 
  select(-variable)


scen1_dalys <- daly_list[["Scenario1.dr3.main"]] %>% as.data.frame() %>% 
  select(param_set, year, dalys) %>% 
  rename(dalys.scen1 = dalys) %>% 
  # Cumulative DALYs, per parameter set
  group_by(param_set) %>% arrange(year) %>% 
  mutate(cum_dalys.scen1 = cumsum(dalys.scen1)) %>% 
  mutate(year = as.integer(year)) # needed for merge to work

scen2_dalys <- daly_list[["Scenario2.dr3.main"]] %>% as.data.frame() %>% 
  select(param_set, year, dalys) %>% 
  rename(dalys.scen2 = dalys) %>% 
  # Cumulative DALYs, per parameter set
  group_by(param_set) %>% arrange(year) %>% 
  mutate(cum_dalys.scen2 = cumsum(dalys.scen2), by = param_set) %>% 
  mutate(year = as.integer(year)) # needed for merge to work

# Calculate NHB over time -- DALYs 

threshold <- 750

NHB_timeDF <- costs_cum_inc_2  %>%  
  left_join(scen1_dalys) %>% 
  left_join(scen2_dalys, by=c("year","param_set")) %>% 
  mutate(inc_dalys = cum_dalys.scen1 - cum_dalys.scen2,
         
         # NHB = dalys averted - incremental costs / threshold
         NHB = inc_dalys - (inc_costs/threshold)) %>% 
  
  group_by(year) %>% 
  summarize(nhb = mean(NHB),
            nhb_min = min(NHB),
            nhb_max = max(NHB))

NHB_varyingDF <- costs_cum_inc_2  %>%  
  left_join(scen1_dalys) %>% 
  left_join(scen2_dalys, by=c("year","param_set")) %>% 
  filter(year==2060) %>% 
  mutate(inc_dalys = cum_dalys.scen1 - cum_dalys.scen2) %>%
  merge(as.data.frame(threshold_varying), .) %>% 
  mutate(NHB = inc_dalys - (inc_costs/threshold_varying)) %>% 
  group_by(threshold_varying) %>% 
  summarize(nhb = mean(NHB),
            nhb_min = min(NHB),
            nhb_max = max(NHB)) 


# -------------------------------------------------------------------------------------------

# Calculate NHB with varying CEA treshold  -- DALYs  (need to also make with QALYs for appendix)

threshold_varying <- 1:1600
CEA_intercept <- 102

# PLOT

Figure3a <- ggplot(data=NHB_varyingDF, aes(x=threshold_varying, y=nhb/1e6)) +   theme_cowplot() +
  
  # Background lines
   geom_hline(yintercept=seq(-10, 10, 1), color="grey90")+
   geom_vline(xintercept=seq(0, 1500, 250), color="grey90")+
  
  # Data   
  geom_ribbon(aes(ymin=nhb_min/1e6, ymax=nhb_max/1e6), alpha = 0.3) +
  geom_line(color = "#00468BFF", size=1) +   # Lancet Colors
  # confidence interval

  
  # Set limits
  xlim(0,1600) +
  ylim(-10, 10) +
  
  # Labels/Titles
  xlab("Cost-Effectiveness Threshold ($)") +
  ylab("Net Health Benefits, in millions \n (DALYs averted)") +
  ggtitle("A) Net Health Benefits by Cost-Effectiveness Threshold, 2060") +
  
  # Marcation lines
  geom_hline(yintercept=0, color="black", linetype="dashed")+
  geom_vline(xintercept=750, color="#AD002AFF", linetype="dashed")+
  geom_vline(xintercept=CEA_intercept, color="#AD002AFF", linetype="dashed")+
  
  # Themes

  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18),
        title = element_text(size = 20)) +
  
  # Add notes
  annotate(geom = "label", x = 970, y = -7, label = "HIV Investment Case \n Threshold = $750", size= 6) +
  annotate(geom = "label", x = 320, y = -7, label = paste0("Incremental Cost per \n DALY averted = $102"), size= 6)

ggsave(Figure3a, file = paste0(cea_path,"figures/Figure3a.pdf"), width = 12, height = 10)

##=================================================================================================

# Net Health Benefits over time

threshold <- 750
time_intercept <- min(NHB_timeDF[NHB_timeDF$nhb>0,"year"]) 

# PLOT

Figure3b <- ggplot(data=NHB_timeDF, aes(x=year, y=nhb/1e6)) +
  
  # Set limits
 scale_x_continuous(limits = c(2020, 2060), breaks = seq(2020, 2060, 10))+
 ylim(-10, 10) +
  
  # Background lines
  geom_hline(yintercept=seq(-10, 10, 1), color="grey90")+
  geom_vline(xintercept=seq(2020,2060, 5), color="grey90")+
  
  
  # Data   
  geom_line(color = "#00468BFF", size=1) +   # Lancet Colors
  # confidence interval
  geom_ribbon(aes(ymin=nhb_min/1e6, ymax=nhb_max/1e6), alpha = 0.3) +
  
  # Labels/Titles
  xlab("Year") +
  ylab("Net Health Benefits, in millions \n (DALYs averted)") +
  ggtitle("B) Net Health Benefits over time, using a threshold of $750") +
  
  # Marcation lines
  geom_hline(yintercept = 0, color="black", linetype="dashed") +
  geom_vline(xintercept = 750, color="#AD002AFF", linetype="dashed") +
  geom_vline(xintercept = time_intercept, color="#AD002AFF", linetype="dashed") +
  
  # Themes
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18),
        title = element_text(size = 20)) +
  
  # Add notes
  annotate(geom = "label", x = 2030, y = -7, label = "First year of positive \n Net Health Benefits = 2023", size= 6) 

ggsave(Figure3a, file = paste0(cea_path,"figures/Figure3b.pdf"), width = 12, height = 8)

# -----------------------------------------------------------------------------------------


fig3 <- ggarrange(Figure3a, Figure3b,  ncol = 1,
                     label.x = 1, label.y = 1)

ggsave(fig3, file = paste0(cea_path,"figures/Figure3.png"), width = 12, height = 10)


