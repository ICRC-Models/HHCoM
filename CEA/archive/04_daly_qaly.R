########################################################################
## Allen Roberts / Mita Sahu
## January 2021 / July 2021
## Calculate DALYs and QALYs from DO ART model ouput scenarios
########################################################################


library(tidyverse)
library(readxl)

## Weights
weights <- expand.grid("hiv_status" = c(0, 1),
                       "cd4_count" = c(0, 1, 2, 3),
                       "art_status" = c(0, 1))
weights <- weights[!(weights$hiv_status == 1 & weights$cd4_count == 0 & weights$art_status == 0), ] ## This combination doesn't exist

## Disability weights
weights$disability_weight[weights$cd4_count == 1] <- 0.582
weights$disability_weight[weights$cd4_count == 2] <- 0.274
weights$disability_weight[weights$cd4_count == 3] <- 0.078
weights$disability_weight[weights$art_status == 1] <- 0.078
weights$disability_weight[weights$hiv_status == 0] <- 0

## Quality of life weights
weights$qaly_weight[weights$cd4_count == 1] <- 0.7
weights$qaly_weight[weights$cd4_count == 2] <- 0.82
weights$qaly_weight[weights$cd4_count == 3] <- 0.94
weights$qaly_weight[weights$art_status == 1] <- 0.94
weights$qaly_weight[weights$hiv_status == 0] <- 1

#================================================================================================================================

# Calculate DALYs and QALYs

daly_list <- vector(mode = "list")
daly_list_summary <- vector(mode = "list")

for (v in names(version)) {  
    
  for (x in scenarios) {
      
      for (d in names(discount_rate)) {
        
        dat <- read_excel(paste0(main_path, x, "/do_art_daly_template_4allen_S", substring(x, first = 9, last = 1e6),".xlsx")) %>% 
          
          filter(param_set!=0) %>% # this is the mean
          
          left_join(weights, by = c("hiv_status", "art_status", "cd4_count")) %>% 
          
          mutate(discount_factor = (1+discount_rate[d])^(year - 2020)) %>%  ## Note that this leaves the first year as undiscounted)
          
          # YLDs
          
          mutate(ylds = (pop_size - num_deaths)*disability_weight/discount_factor) %>% 
        
          # YLLs
          
          mutate(end_year = pmin(hrzn, year + (80-age)))
        
          # Discount YLLs
        
          for(i in 1:nrow(dat)) {
            dat$ylls[i] <- dat$num_deaths[i]*sum((1/dat$discount_factor[i]) * (1/(1+discount_rate))^(0:(dat$end_year[i] - dat$year[i])))
          }
        
        ## QALYs
        
        dat$qalys <- (dat$pop_size-dat$num_deaths)*dat$qaly_weight/dat$discount_factor
        
        ## Summarize
        daly <- dat %>%
          group_by(param_set) %>% 
          group_by(year) %>% 
          summarize(ylls = sum(ylls),
                    ylds = sum(ylds),
                    dalys = sum(ylds) + sum(ylls),
                    qalys = sum(qalys))
        
        temp <- dat %>%
          group_by(param_set) %>%
          summarize(ylls = sum(ylls),
                    ylds = sum(ylds),
                    dalys = sum(ylds) + sum(ylls),
                    qalys = sum(qalys))
        
        summary <- data.frame(param = c("ylls", "ylds", "dalys", "qalys"),
                             mean = c(mean(temp$ylls), mean(temp$ylds), mean(temp$dalys), mean(temp$qalys)),
                             min =  c(min(temp$ylls), min(temp$ylds), min(temp$dalys), min(temp$qalys)),
                             max = c(max(temp$ylls), max(temp$ylds), max(temp$dalys), max(temp$qalys)))
        
        daly_list[[length(daly_list) + 1]] <- daly
        daly_list_summary[[length(daly_list_summary) + 1]] <- summary
        
        names(daly_list)[length(daly_list)] <- paste(x, d, v, sep = ".")
        names(daly_list_summary)[length(daly_list_summary)] <- paste(x, d, v, sep = ".")
        
    }
  }
}

