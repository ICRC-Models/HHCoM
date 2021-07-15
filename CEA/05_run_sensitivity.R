# Sensitivity Analysis
# MSahu
# July 14, 2021

# RUN 

icerMAT <- matrix(NA, nrow=nrow(snDF), ncol=3)

for (j in 1:nrow(snDF)) {
  
  snDF[j,"ON_or_OFF"] = T
  
  print(snDF[j, "sn_name_full"])
  
  source(paste0(cea_path, "02_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
  source(paste0(cea_path, "04_icer_bia.R"))
  
  snDF[j,"ON_or_OFF"] = F
  
  icerMAT[j, ] <- icer
}

snDF$icer_mean <- icerMAT[,1]
snDF$icer_min <- icerMAT[,2]
snDF$icer_max <- icerMAT[,3]

tornadoDF <- snDF %>% 
  select(sn_name, bound_abb, icer_mean, icer_min, icer_max) %>% 
  mutate(icer = ifelse(bound_abb == "lb",icer_min, icer_max)) %>% 
  reshape2::dcast(sn_name ~ bound_abb, value.var = "icer") %>% 
  mutate(level = -row_number())
#  mutate(level = floor(-row_number()/2))


# original value of output
base.value <- icer[1,1]

# get order of parameters according to size of intervals
# (I use this to define the ordering of the factors which I then use to define the positions in the plot)
order.parameters <- tornadoDF %>% 
  mutate(Parameter=factor(x=sn_name, levels=sn_name)) %>%
  select(Parameter) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.95

# get data frame in shape for ggplot and geom_rect
df.2 <- df %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value', Lower_Bound:Upper_Bound) %>%
  # just reordering columns
  select(Parameter, type, output.value, UL_Difference) %>%
  # create the columns for geom_rect
  mutate(Parameter=factor(Parameter, levels=order.parameters),
         ymin=pmin(output.value, base.value),
         ymax=pmax(output.value, base.value),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2)


# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
library(ggplot2)
ggplot(data = tornadoDF, 
       aes(ymax=level, ymin=level-1, xmax=ub, xmin=lb, fill = sn_name)) + 
  geom_rect() +
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_vline(xintercept = base.value) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  coord_flip()

