########################################################################
## Allen Roberts
## January 2021
## Calculate DALYs and QALYs from DO ART model ouput scenarios
########################################################################
rm(list = ls())

library(tidyverse)
library(readxl)

## Parameters
discount_rate <- 0.03 ## Set to 0 for undiscounted
time_horizon <- 2060 ## Right now this is inclusive - eg, start of 2020 through end of 2060. This is actually 41 years, not 40 years
time_step <- 1 ## years

## Load data and format
s1 <- read_excel("do_art_daly_template_4allen_S1.xlsx")
s2 <- read_excel("do_art_daly_template_4allen_S2.xlsx")
s3 <- read_excel("do_art_daly_template_4allen_S3.xlsx")

s1$scenario <- "Scenario 1"
s2$scenario <- "Scenario 2"
s3$scenario <- "Scenario 3"

dat <- rbind(s1, s2, s3)

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

dat <- dat %>%
  left_join(weights, by = c("hiv_status", "art_status", "cd4_count"))

## Discount factor
dat$discount_factor <- (1+discount_rate)^(dat$year - 2020) ## Note that this leaves the first year as undiscounted

## YLDs
dat$ylds <- (dat$pop_size - dat$num_deaths)*dat$disability_weight*time_step/dat$discount_factor

## YLLs
dat$end_year <- pmin(time_horizon, dat$year + (80-dat$age))
dat$ylls <- NA

for(ii in 1:nrow(dat)) {
  dat$ylls[ii] <- dat$num_deaths[ii]*sum((1/dat$discount_factor[ii])*(1/(1+discount_rate))^(0:(dat$end_year[ii] - dat$year[ii])))
}

## QALYs
dat$qalys <- (dat$pop_size-dat$num_deaths)*dat$qaly_weight*time_step/dat$discount_factor

## Summarize
output <- dat %>%
  group_by(scenario) %>%
  summarize(ylls = sum(ylls),
            ylds = sum(ylds),
            dalys = sum(ylds) + sum(ylls),
            qalys = sum(qalys))

output$ylds_diff <- output$ylds - output$ylds[output$scenario == "Scenario 1"]
output$ylls_diff <- output$ylls - output$ylls[output$scenario == "Scenario 1"]
output$dalys_diff <- output$dalys - output$dalys[output$scenario == "Scenario 1"]
output$qalys_diff <- output$qalys - output$qalys[output$scenario == "Scenario 1"]

## Save
output$discount_rate <- discount_rate
write.csv(output, file = paste0("output_dr", as.integer(100*discount_rate), ".csv"), row.names = FALSE)
