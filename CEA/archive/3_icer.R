# ICER

#MSahu
#July 28, 2020

library(ggplot2)

icer <- mortality %>% rename(mortality_2_1_diff = diff_2_1_mean) %>% select(year, mortality_2_1_diff) %>% 
  left_join(incidence, by ="year") %>%  rename(incidence_2_1_diff = diff_2_1_mean) %>% 
  select(year, incidence_2_1_diff, mortality_2_1_diff) %>% 
  left_join(cost_diff) %>% 
  filter(year>=2016) %>% 
  mutate(cost_death_averted= cost_diff_1_2/(mortality_2_1_diff*100000)) %>% 
  mutate(cost_case_averted=cost_diff_1_2/(incidence_2_1_diff*100))
  
ggplot(data=icer, aes(x=year, y=cost_death_averted)) +
  geom_line() +
  xlab("Year") +
  ylab("Cost per death averted") +
  theme_bw()

ggplot(data=icer, aes(x=year, y=cost_case_averted)) +
  geom_line() +
  xlab("Year") +
  ylab("Cost per case averted") +
  theme_bw()