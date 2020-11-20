
# PLOTS

# Annual dalys (total)

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_scen1.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual DALYs for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,350)

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_scen2.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual DALYs for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,350)

# Cumulative dalys (total)

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_cum_scen1.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/100000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative DALYs for KZN / 100k") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_cum_scen2.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/100000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative DALYs for KZN / 100k") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) 
