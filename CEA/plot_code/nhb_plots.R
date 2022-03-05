# ICER and  CEA planes ; time plots
# September 10, 2020
# MSahu

library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(reshape2)

##=================================================================================================

# Net Health Benefits over time

threshold <- 750

costs_cum_inc_2 <-  costs_cum_incr_list$Scenario2.total.dr3.main %>% 
  select(year,mean, min, max) 


dalys_cum_averted_2 <- read.csv(paste0(cea_path,"effects/daly/daly_cum_averted_raw_scen2.csv")) %>% 
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable)

# Calculate NHB so positive is good

NHB_daly <- costs_cum_inc_2  %>%  left_join(dalys_cum_averted_2, by=c("year","set_name")) %>% 
  mutate(nhb = (dalys - (costs/threshold))) %>% 
  group_by(year) %>% 
  summarise_at(vars(nhb), funs(mean, min, max)) %>% 
  rename(nhb = mean,
         nhb_min = min,
         nhb_max = max)


# PLOT

# "Intercept" is defined BELOW

Figure2a <- ggplot(data=NHB_daly, aes(x=year, y=nhb/1000)) +
  
  # Set limits
  xlim(2020, 2060) +
  ylim(-15000, 15000) +
  
  # Background lines
  geom_hline(yintercept=seq(-15000, 15000, 2500), color="grey90")+
  geom_vline(xintercept=seq(2020, 2060, 10), color="grey90")+
  
  # Data   
  geom_line(color = "#00468BFF", size=1) +   # Lancet Colors
  # confidence interval
  geom_ribbon(aes(ymin=nhb_min/1000, ymax=nhb_max/1000), alpha = 0.3) +
  
  # Labels/Titles
  xlab("Year") +
  ylab("Net Health Benefits, in thousands \n (Health gains measured using DALYs averted)") +
  ggtitle("Net Health Benefits for Community-Based ART over time, using a threshold of 750 USD") +
  
  # Marcation lines
  geom_hline(yintercept=0, color="black", linetype="dashed") +
  geom_vline(xintercept=750, color="#AD002AFF", linetype="dashed") +
  geom_vline(xintercept=intercept, color="#AD002AFF", linetype="dashed") +
  
  # Themes
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

ggsave(Figure2a, file = paste0(cea_path,"figures/Figure2a.pdf"), width = 12, height = 8)




###############################################################################################################

# NHB with varying CEA treshold  -- DALYs  (need to also make with QALYs for appendix)


threshold <- 1:1600

costs_cum_inc_2 <-  costs_cum_incr_list$Scenario2.total.dr3.main %>% 
  select(year,mean, min, max) 


dalys_cum_averted_2 <- read.csv(paste0(cea_path,"effects/daly/daly_cum_averted_raw_scen2.csv")) %>% 
  filter(year==2060) %>% 
  select(-c(mean, min, max)) %>% 
  melt(id="year",value.name="dalys_averted") %>% 
  rename(set = variable)

merged <- costs_cum_inc_2 %>%  left_join(dalys_cum_averted_2, by=c("year","set")) 


NHB_daly <- merge(as.data.frame(threshold), merged) %>% 
  mutate(nhb = (dalys_averted - (costs/threshold))) %>% 
  group_by(year, threshold) %>% 
  summarise_at(vars(nhb), funs(mean, min, max)) %>% 
  rename(nhb = mean,
         nhb_min = min,
         nhb_max = max)

intercept <- pull(NHB_daly[abs(NHB_daly$nhb)==min(abs(NHB_daly$nhb)),"threshold"]) 


# PLOT

Figure2b <- ggplot(data=NHB_daly, aes(x=threshold, y=nhb/1000)) +
  
  # Set limits
  xlim(55,1600) +
  ylim(-15000, 15000) +
  
  # Background lines
  geom_hline(yintercept=seq(-15000, 15000, 2500), color="grey90")+
  geom_vline(xintercept=seq(0, 1500, 250), color="grey90")+

  # Data   
  geom_line(color = "#00468BFF", size=1) +   # Lancet Colors
  # confidence interval
  geom_ribbon(aes(ymin=nhb_min/1000, ymax=nhb_max/1000), alpha = 0.3) +
  
  # Labels/Titles
  xlab("Cost-Effectiveness Threshold (in USD)") +
  ylab("Net Health Benefits, in thousands \n (Health gains measured using DALYs averted)") +
  ggtitle("Net Health Benefits for Community-Based ART, 2060") +
  
  # Marcation lines
  geom_hline(yintercept=0, color="black", linetype="dashed")+
  geom_vline(xintercept=750, color="#AD002AFF", linetype="dashed")+
  geom_vline(xintercept=intercept, color="#AD002AFF", linetype="dashed")+

  # Themes
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  
  # Add notes
  annotate(geom = "label", x = 950, y = -10000, label = "HIV Investment Case \n Threshold = $750", size= 6) +
  annotate(geom = "label", x = 300, y = -10000, label = paste0("Incremental Cost per \n DALY averted = $148"), size= 6)
  
ggsave(Figure2b, file = paste0(cea_path,"figures/Figure2b.pdf"), width = 12, height = 8)