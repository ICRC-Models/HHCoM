# QALY plots

# Annual qalys (total)

qalys <- read.csv(paste0(cea_path,"effects/qaly/qaly_scen1.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys, aes(x=year, y=qalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual QALYs for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,800)

qalys <- read.csv(paste0(cea_path,"effects/qaly/qaly_scen2.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys, aes(x=year, y=qalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual QALYs for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,800)

# Cumulative qalys (total)

qalys <- read.csv(paste0(cea_path,"effects/qaly/qaly_cum_scen1.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys, aes(x=year, y=qalys/100000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative QALYs for KZN / 100k") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,2500)

qalys <- read.csv(paste0(cea_path,"effects/qaly/qaly_cum_scen2.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="qalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=qalys, aes(x=year, y=qalys/100000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Cumulative QALYs for KZN / 100k") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,2500)

