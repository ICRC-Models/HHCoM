# Calculate DALYs
# MSahu
# October 6, 2020

# DALY Calculation: YLL + (YLD * Disability Weight)
# This script calculates average annual DALYs averted between interventions


# Note that health benefits for INCIDENT cases are not included. - ie. 
# assumption is that health benefits show up the next year.

# Note that this is a Ctrl - F find and replace QALY -- > DALY from script #3

##############################################################################


# Set up utility weight parameters

wt_hiv_neg <- 0
wt_cd4_350 <- 0.012
wt_cd4_250_350 <- 0.274
wt_cd4_200 <- 0.582
wt_art <- 0.078
wt_dead <- 1

# Current strategy - calculate on an annual basis ; from the START of the year and don't include incident cases.
#1. Number of prevalent cases * utility weight (by CD4+ count) 
#2. Number of HIV negative cases * utility weight  [can be calculated as 1 - prevalence]
#3. Sum utilities
#4. Discount

##############################################################################

# IMPORT POPULATION (males and females combined, age 15-79)

##############################################################################

for (x in 1:3) {
  
  pop <- read.csv(paste0(main_path,"Scenario",x,"/PopulationSize_combined_aged15-79.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("pop_scen",x),pop)
  
}
rm(pop)

##############################################################################

# IMPORT PREVALENCE (males and females combined), by CD4+ count and ART status

##############################################################################

# Prevalence on ART

for (x in 1:3) {
  
  prev_ART<- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_onART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_ART_scen",x),prev_ART)
  
}

# Prevalence CD4 below 200, not on ART

for (x in 1:3) {
  
  prev_200 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_below200_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_200_scen",x),prev_200)
  
}

rm(lists=ls(pattern="mort")

# Prevalence CD4 200-350, not on ART

for (x in 1:3) {
  
  prev_200_350 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_200-350_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_200_350_scen",x),prev_200_350)
  
}

# Prevalence CD4 350+, not on ART 

for (x in 1:3) {
  
  prev_350 <- read.csv(paste0(main_path,"Scenario",x,"/HIV_prevalence_combined_aged15-79_CD4_350plus_noART.csv"), header=F) %>% 
    setNames(paste0("s",-3:25)) %>% 
    dplyr::rename(year=1,
                  mean=2,
                  min=3,
                  max=4) %>% 
    filter(year>=2020)
  
  assign(paste0("prev_350_scen",x),prev_350)
  
}


# MULTIPLY BY UTILITY WEIGHTS 

for (x in 1:3) {
  
  hiv_pos <- (get(paste0("prev_ART_scen",x))[,-1] * wt_art +    
                get(paste0("prev_200_scen",x))[,-1] * wt_cd4_200 +
                get(paste0("prev_200_350_scen",x))[,-1] * wt_cd4_250_350 +
                get(paste0("prev_350_scen",x))[,-1] * wt_cd4_350  ) %>% 
    mutate(year=2020:2060) %>% 
    select(year, everything()) 
  
  assign(paste0("hiv_pos_scen",x),hiv_pos)
  
}

##############################################################################

# DEATHS (DISABILITY WEIGHT = 1)

##############################################################################

#MULTIPLY BY DISABILITY WEIGHT

for (x in 1:3) {
  
  mortality <- ((( read.csv(paste0(main_path,"Scenario",x,"/HIV_mortality_combined_aged15-79.csv"), header=F)[96:136,-1])    # import 2020-2060 only
                          * wt_dead )      # disability weight for death
                          / 100000 ) %>%   # mortality rates are given per 100k
    setNames(paste0("s",-2:25)) %>% 
    dplyr::rename(mean=1,
                  min=2,
                  max=3) %>% 
    mutate(year=2020:2060) %>% 
    select(year, everything())
  
  
  assign(paste0("mortality",x),mortality)
  
}


# Remove excess DFs

rm(list=ls(pattern="prev"))


##############################################################################

# Annual DALYs - sum HIV positive and HIV negative, then discount

##############################################################################


for (x in 1:3) {
  
  daly <- ( ( get(paste0("mortality",x))[,-1] + get(paste0("hiv_pos_scen",x))[,-1] )
            
            # MULTIPLY BY POPULATION
            * get(paste0("pop_scen",x))[,-1] ) %>% 
    
    # DISCOUNT
    mutate(year_discount=0:40,   # set 2020 to Year 0
           discount_amt=discount(year_discount=year_discount,
                                 discount_rate=discount_rate)) %>% 
    mutate_at(1:28,~discounter(.,discount_amt)) %>% 
    select(-year_discount,-discount_amt) %>% 
    
    mutate(year=2020:2060) %>% 
    select(year,everything()) 
  
  assign(paste0("daly_scen",x),daly)
  
}

rm(daly)


# Cumulative 

for (x in 1:3) {
  
  daly <- get(paste0("daly_scen",x))[,-1] %>% 
            
    transmute_at(1:28, ~cumsum(.)) %>% 
    
    mutate(year=2020:2060) %>% 
    select(year,everything()) 
  
  assign(paste0("daly_cum_scen",x),daly)
  
}

rm(daly)

###################################################################################################################

# Annual raw dalys averted

daly_averted_raw_scen2 <- (daly_scen1[,-1] - daly_scen2[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

daly_averted_raw_scen3 <- (daly_scen1[,-1] - daly_scen3[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Annual percent dalys averted

daly_averted_pct_scen2 <- (((daly_scen1[,-1] - daly_scen2[,-1] )/daly_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

daly_averted_pct_scen3 <- (((daly_scen1[,-1] - daly_scen3[,-1] )/daly_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative raw dalys averted

daly_cum_averted_raw_scen2 <- ( daly_cum_scen1[,-1] - daly_cum_scen2[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

daly_cum_averted_raw_scen3 <- ( daly_cum_scen1[,-1] - daly_cum_scen3[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative percent dalys averted

daly_cum_averted_pct_scen2 <- (((daly_cum_scen1[,-1] - daly_cum_scen2[,-1] )/daly_cum_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

daly_cum_averted_pct_scen3 <- (((daly_cum_scen1[,-1] - daly_cum_scen3[,-1] )/daly_cum_scen1[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

###################################################################################################################

# Export CSVs

csv_list <- list("daly_scen1",
                 "daly_scen2",
                 "daly_scen3",
                 "daly_cum_scen1",
                 "daly_cum_scen2",
                 "daly_cum_scen3",
                 "daly_averted_raw_scen2",
                 "daly_averted_raw_scen3",
                 "daly_averted_pct_scen2",
                 "daly_averted_pct_scen3",
                 "daly_cum_averted_raw_scen2",
                 "daly_cum_averted_raw_scen3",
                 "daly_cum_averted_pct_scen2",
                 "daly_cum_averted_pct_scen3")


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"effects/daly/",x,".csv"), row.names = F))

# Remove excess DFs

rm(list=ls())

#############################################################################################

# PLOTS



# Annual DALYs (total)

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_scen1.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual dalys for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,40)

dalys <- read.csv(paste0(cea_path,"effects/daly/daly_scen2.csv")) %>% 
  # reshape long
  reshape2::melt(id="year",value.name="dalys") %>% 
  rename(set_name=variable) %>% 
  mutate(size=ifelse(set_name %in% c("mean","min","max"),"bold","normal"))

ggplot(data=dalys, aes(x=year, y=dalys/10000,group=set_name)) +
  geom_line(aes(color=set_name,size=size)) +
  scale_size_manual(values=c(1.5,.1)) +
  xlab("Year") +
  ylab("Annual dalys for KZN / 10000") +
  theme_cowplot() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=18)) +
  ylim(0,40)

# Cumulative DALYs (total)

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
        axis.text=element_text(size=18)) +
  ylim(0,70)

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
        axis.text=element_text(size=18)) +
  ylim(0,70)

#  Scenario 2


costs <- read.csv(paste0(cea_path,"costs/inc_costs_cum_pc_scen2.csv")) %>% 
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