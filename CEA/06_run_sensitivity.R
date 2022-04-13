# Sensitivity Analysis
# MSahu
# July 14, 2021

# -----------------------------------------------------------------------------------------------------

library(ggplot2)

# First run main script 

source(paste0(cea_path, "helper.R"))
source(paste0(cea_path, "01_setup_sensitivity.R"))
source(paste0(cea_path, "02_cases_deaths.R"))
source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

# original ICER mean
icer_base <- mean(icer)


# RUN SENSITIVITY

icerMAT <- matrix(NA, nrow=nrow(snDF), ncol=3)

for (j in 1:nrow(snDF)) {
  
  snDF[j,"ON_or_OFF"] = T
  
  print(snDF[j, c("sn_name_full", "bound")])
  
  source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
  source(paste0(cea_path, "05_icer_bia.R"))
  
  snDF[j,"ON_or_OFF"] = F
  
  icerMAT[j, ] <- c(mean(icer), min(icer), max(icer))
}

snDF$icer_mean <- icerMAT[,1]
snDF$icer_min <- icerMAT[,2]
snDF$icer_max <- icerMAT[,3]

# -----------------------------------------------------------------------------------

# TWO-WAY sensitivity 

# UPPER: upper bound for cb_ART; lower bound for clinic ART

snDF[6,"ON_or_OFF"] = T
snDF[7,"ON_or_OFF"] = T

print(snDF[6, c("sn_name_full", "bound")])
print(snDF[7, c("sn_name_full", "bound")])

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

upper <- icer

snDF[6,"ON_or_OFF"] = F
snDF[7,"ON_or_OFF"] = F

# LOWER: lower bound for cb_ART; upper bound for clinic ART


snDF[5,"ON_or_OFF"] = T
snDF[8,"ON_or_OFF"] = T

print(snDF[5, c("sn_name_full", "bound")])
print(snDF[8, c("sn_name_full", "bound")])

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

lower <- icer

snDF[5,"ON_or_OFF"] = F
snDF[8,"ON_or_OFF"] = F

# Create lists to bind with dataframe for Tornado Diagram
two_way <- list("Two-way: Cost of Community and Clinic ART", mean(lower), mean(upper), min(lower), min(upper), max(lower), max(upper))
two_way.lb <- list("Two-way: Cost of Community and Clinic ART", "lb", mean(lower), min(lower), max(lower))
two_way.ub <- list("Two-way: Cost of Community and Clinic ART", "ub", mean(upper), min(upper), max(upper))

# --------------------------------------------------------------------------------------

# Percentage of people on cb-ART 
# This impacts the scalars for Scenario 2 -- lines 425-447 of 03_costs Script
# Parameters for the bounds are stored in: paste0(cea_path,"parameters/art_parameters.xls"), sheet = "Table_3")

# LOWER

sn_pct_cbART_main = F
sn_pct_cbART_lo = T

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

lower <- icer

# UPPER

sn_pct_cbART_lo = F
sn_pct_cbART_hi = T

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

upper <- icer

# Create lists to bind with dataframe for Tornado Diagram

name_cbART_rate <- "Proportion at Community ART Cost (0-100%)"
cbART_rate <- list(name_cbART_rate, mean(lower), mean(upper), min(lower), min(upper), max(lower), max(upper))
cbART_rate.lb <- list(name_cbART_rate, "lb", mean(lower), min(lower), max(lower))
cbART_rate.ub <- list(name_cbART_rate, "ub", mean(upper), min(upper), max(upper))

# Clean Up

sn_pct_cbART_main = T
sn_pct_cbART_hi = F

# -----------------------------------------------------------------------------------

# TWO-WAY sensitivity  (with cb-ART rate)

# UPPER: upper bound for cb_ART; upper bound for cb-ART factor

 snDF[6,"ON_or_OFF"] = T
 sn_pct_cbART_main = F
 sn_pct_cbART_hi = T

 source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
 source(paste0(cea_path, "05_icer_bia.R"))

 upper <- icer

 snDF[6,"ON_or_OFF"] = F


# LOWER: original

sn_pct_cbART_hi = F
sn_pct_cbART_main = T

snDF[5,"ON_or_OFF"] = T

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

lower <- icer

snDF[5,"ON_or_OFF"] = F

# Create lists to bind with dataframe for Tornado Diagram
two_way2 <- list("Two-way: Community ART Cost and Proportion (0-100%)", mean(lower), mean(upper), min(lower), min(upper), max(lower), max(upper))
two_way2.lb <- list("Two-way: Community ART Cost and Proportion (0-100%)", "lb", mean(lower), min(lower), max(lower))
two_way2.ub <- list("Two-way: Community ART Cost and Proportion (0-100%)", "ub", mean(upper), min(upper), max(upper))


#-----------------------------------------------------------------------------------------

# Discount rate (in 05_icer_bia Script)

# LOWER

sn_dr_main = F
sn_dr_lo = T

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

lower <- icer

# UPPER
sn_dr_lo = F
sn_dr_hi = T

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

upper <- icer

# Create lists to bind with dataframe for Tornado Diagram
name_disc_rate <- paste0("Annual discount rate (",discount_rate[1]*100,"-", discount_rate[3]*100,"%)")
disc_rate <- list(name_disc_rate, mean(lower), mean(upper), min(lower), min(upper), max(lower), max(upper))
disc_rate.lb <- list(name_disc_rate, "lb", mean(lower), min(lower), max(lower))
disc_rate.ub <- list(name_disc_rate, "ub", mean(upper), min(upper), max(lower))

# Clean Up
sn_dr_hi = F
sn_dr_main = T

# ---------------------------------------------------------------------------------------

# TORNADO PLOT

parametersDF <- snDF %>% 
  select(sn_name_full, bound_abb, icer_mean, icer_min, icer_max) %>% setDT() %>% 
  data.table::dcast(sn_name_full ~ bound_abb, value.var = c("icer_mean", "icer_min", "icer_max")) %>%
  rbind(two_way, 
        #two_way2, cbART_rate, 
        disc_rate) %>% 
  mutate(UI_min = apply(.[,2:7], 1, FUN = min),
         UI_max = apply(.[,2:7], 1, FUN = max)) %>% 
  mutate(UI_diff = abs(icer_mean_ub - icer_mean_lb))

# get order of parameters according to size of intervals
# (I use this to define the ordering of the factors which I then use to define the positions in the plot)
order.parameters <- parametersDF %>% arrange(-UI_diff) %>% 
  mutate(Parameter=factor(x=sn_name_full, levels=sn_name_full)) %>%
  select(Parameter) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.8

# get data frame in shape for ggplot and geom_rect
tornadoDF <- snDF %>% 
  select(sn_name_full, bound_abb, icer_mean, icer_min, icer_max) %>% 
  rbind(two_way.lb, two_way.ub, disc_rate.lb, disc_rate.ub) %>% 
  left_join(parametersDF, by = "sn_name_full") %>% 
  select(sn_name_full, bound_abb, icer_mean, icer_min, icer_max, UI_min, UI_max) %>%
  mutate(icer.base = icer_base) %>% 
  # create the columns for geom_rect
  mutate(Parameter=factor(sn_name_full, levels=order.parameters),
         ymin=pmin(icer_mean, icer.base),
         ymax=pmax(icer_mean, icer.base),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2) 

tornadoDF_appendix <- snDF %>% 
  select(sn_name_full, bound_abb, icer_mean, icer_min, icer_max) %>% 
  rbind(two_way.lb, two_way.ub, 
        two_way2.lb, two_way2.ub, 
        cbART_rate.lb, cbART_rate.ub, disc_rate.lb, disc_rate.ub) %>% 
  left_join(parametersDF, by = "sn_name_full") %>% 
  select(sn_name_full, bound_abb, icer_mean, icer_min, icer_max, UI_min, UI_max) %>%
  mutate(icer.base = icer_base) %>% 
  # create the columns for geom_rect
  mutate(Parameter=factor(sn_name_full, levels=order.parameters),
         ymin=pmin(icer_mean, icer.base),
         ymax=pmax(icer_mean, icer.base),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2) 


# PLOT
# (use scale_x_continuous to change labels in y axis to name of parameters)

darjeeling_orange = "#F98400"
darjeeling_yellow = "#F2AD00"
darjeeling_green = "#00A08A"
darjeeling_blue = "#5BBCD6"
darjeeling_red = "#FF0000"
lancet_red = "#AD002AFF"
lancet_green = "#42B540FF"
lancet_blue = "#00468BFF"

if (hrzn %in% c(2030, 2060)) {tornado_UB = 800}
if (hrzn == 2045) {}

tornado_UB = 800
tornado_LB = -100

tornado <- ggplot(data = tornadoDF) + 
  geom_hline(yintercept = seq(tornado_LB,tornado_UB, by = 100), size = 1, color = "gray90") +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = bound_abb)) +
  geom_errorbar(aes(ymin = UI_min, ymax = UI_max, x = as.numeric(Parameter)), width = 0.2) +
  theme_classic() + ylab("ICER per DALY averted ($)") +
  scale_fill_manual(values = c(darjeeling_blue, lancet_red), labels = c("Lower bound", "Upper bound")) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.title.x = element_text(size = 30), 
        axis.text = element_text(size=30),
        legend.text = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_hline(yintercept = threshold, linetype = "dashed", size = 2, color = lancet_blue) +
  geom_hline(yintercept = 0, size = 1) +
  geom_hline(yintercept = icer_base) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  scale_y_continuous(limits = c(tornado_LB, tornado_UB),
                  breaks = seq(tornado_LB, tornado_UB, by = 100),
                  labels = seq(tornado_LB, tornado_UB, by = 100))+
  annotate(geom = "text", x = 4, y = 550, label = paste0("Cost-Effectiveness Threshold \n = $750 per DALY averted"), 
           size= 10, color = lancet_blue) + 
  coord_flip()

# Save

ggsave(plot = tornado, file = paste0(dir, "CEA/figures/Figure4_",hrzn,".pdf"), device = "pdf",
       width = 25, height = 10)
ggsave(plot = tornado, file = paste0(dir, "CEA/figures/Figure4_",hrzn,".eps"), device = "eps",
       width = 25, height = 10)
