# Set up of parameters ------------------------------------------

fixedpars = list()

fixedpars$sickpop = 1250

fixedpars$pr_seek = 0.8

fixedpars$sick = 1000 

fixedpars$pr_comp = c(1, 0.9, 1) # reduces hospitalization by 10% (hence reducing mortality by 10%) or pay twice as much to reduce only mortality by 20%
# make the sensitivity analysis show that it depends a lot on hospitalization
fixedpars$pr_simple = 1-fixedpars$pr_comp
fixedpars$pr_recover = c(0.5, 0.5, 0.6)
fixedpars$pr_death = 1-fixedpars$pr_recover

fixedpars$cost_rx = c(0, 300, 600)
fixedpars$cost_hosp = 500

fixedpars$yld_comp = 0.1
fixedpars$yld_simple = 0.05

fixedpars$yll = 10
