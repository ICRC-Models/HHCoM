# Calculate QALYs
# MSahu
# September 9, 2020

# QALY Calculation: Years of Life * Utility Value
# This script calculates average annual QALYs gained between interventions


# TO DO:
# Check if CD4+ and ART status of incident cases for model is available


##############################################################################


# Set up utility weight parameters

wt_hiv_neg <- 1
wt_cd4_350 <- 0.94
wt_cd4_250_350 <- 0.82
wt_cd4_200 <- 0.7
wt_art <- 0.94
wt_dead <- 0

# Current strategy - calculate on an annual basis ; from the START of the year and don't include incident cases.
#1. Number of prevalent cases * utility weight (by CD4+ count) 
#2. Number of HIV negative cases * utility weight  [can be calculated as 1 - prevalence]
#3. Sum utilities
#4. Discount

# Then Ask Paul Revill - because actually should be cumulative ! (or is there a way to track births?)


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

# HIV-NEGATIVE (males and females combined), by CD4+ count and ART status

##############################################################################

#MULTIPLY BY UTILITY WEIGHT

for (x in 1:3) {
  
  hiv_neg <- ((1- (get(paste0("prev_ART_scen",x))[,-1] +    # HIV-negative is 1 - HIV Prevalence
                get(paste0("prev_200_scen",x))[,-1] +
                get(paste0("prev_200_350_scen",x))[,-1] +
                get(paste0("prev_350_scen",x))[,-1] )) 
                      * wt_hiv_neg )  %>%   #multiply by utility weight (in this case = 1)
  mutate(year=2020:2060) %>% 
  select(year, everything()) 
  
  assign(paste0("hiv_neg_scen",x),hiv_neg)
  
}

# Remove excess DFs

rm(list=ls(pattern="prev"))
rm(hiv_neg, hiv_pos)

##############################################################################

# Annual QALYs - sum HIV positive and HIV negative, then discount

##############################################################################


for (x in 1:3) {
  
  qaly <- ( ( get(paste0("hiv_neg_scen",x))[,-1] + get(paste0("hiv_pos_scen",x))[,-1] )
  
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
  
  assign(paste0("qaly_scen",x),qaly)
  
}

rm(qaly)


# Cumulative 

for (x in 1:3) {
  
  qaly <- ( ( get(paste0("hiv_neg_scen",x))[,-1] + get(paste0("hiv_pos_scen",x))[,-1] )
            
            # MULTIPLY BY POPULATION
            * get(paste0("pop_scen",x))[,-1] )  %>% 
    
    # DISCOUNT
    mutate(year_discount=0:40,   # set 2020 to Year 0
           discount_amt=discount(year_discount=year_discount,
                                 discount_rate=discount_rate)) %>% 
    mutate_at(1:28,~discounter(.,discount_amt)) %>% 
    select(-year_discount,-discount_amt) %>% 
    transmute_at(1:28, ~cumsum(.)) %>% 
    
    mutate(year=2020:2060) %>% 
    select(year,everything()) 
  
  assign(paste0("qaly_cum_scen",x),qaly)
  
}


###################################################################################################################

# Annual raw qalys gained

qaly_gained_raw_scen2 <- (qaly_scen2[,-1] - qaly_scen1[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

qaly_gained_raw_scen3 <- (qaly_scen3[,-1] - qaly_scen1[,-1] )  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Annual percent qalys gained

qaly_gained_pct_scen2 <- (((qaly_scen2[,-1] - qaly_scen1[,-1] )/qaly_scen2[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

qaly_gained_pct_scen3 <- (((qaly_scen3[,-1] - qaly_scen1[,-1] )/qaly_scen3[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative raw qalys gained

qaly_cum_gained_raw_scen2 <- ( qaly_cum_scen2[,-1] - qaly_cum_scen1[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

qaly_cum_gained_raw_scen3 <- ( qaly_cum_scen3[,-1] - qaly_cum_scen1[,-1] ) %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

# Cumulative percent qalys gained

qaly_cum_gained_pct_scen2 <- (((qaly_cum_scen2[,-1] - qaly_cum_scen1[,-1] )/qaly_cum_scen2[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

qaly_cum_gained_pct_scen3 <- (((qaly_cum_scen3[,-1] - qaly_cum_scen1[,-1] )/qaly_cum_scen3[,-1] )*100)  %>% 
  mutate(year=2020:2060) %>% 
  select(year, everything()) 

###################################################################################################################

# Export CSVs

csv_list <- list("qaly_scen1",
                 "qaly_scen2",
                 "qaly_scen3",
                 "qaly_cum_scen1",
                 "qaly_cum_scen2",
                 "qaly_cum_scen3",
                 "qaly_gained_raw_scen2",
                 "qaly_gained_raw_scen3",
                 "qaly_gained_pct_scen2",
                 "qaly_gained_pct_scen3",
                 "qaly_cum_gained_raw_scen2",
                 "qaly_cum_gained_raw_scen3",
                 "qaly_cum_gained_pct_scen2",
                 "qaly_cum_gained_pct_scen3")


lapply(csv_list, function(x) write.csv(get(x), file=paste0(cea_path,"effects/qalys/",x,".csv")))

# Remove excess DFs

rm(list=ls())

#############################################################################################

# PLOTS



# Annual qalys (total)

qalys <- read.csv(paste0(cea_path,"effects/qalys/qaly_scen1.csv")) %>% 
  select(-X) %>% 
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

qalys <- read.csv(paste0(cea_path,"effects/qalys/qaly_scen2.csv")) %>% 
  select(-X) %>% 
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

qalys <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_scen1.csv")) %>% 
  select(-X) %>% 
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

qalys <- read.csv(paste0(cea_path,"effects/qalys/qaly_cum_scen2.csv")) %>% 
  select(-X) %>% 
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

