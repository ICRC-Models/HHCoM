
## Parameters
n <- 1000 ## Population size
p <- 0.3 ## HIV prevalence at start of simulation
inc <- 0.02 ## HIV incidence : 2 per 100 PY (approximately 2% per year)
prop_diagnosed <- 0.85 ## First 90 - this is what you would take from SABSSM-V
n_hiv_diag <- n*p*prop_diagnosed ## Number of PLHIV who know their status
n_hiv_undiag <- n*p*(1-prop_diagnosed) ## Number of PLHIV who do not know their status
n_hiv <- n_hiv_diag + n_hiv_undiag ## Total number of PLHIV in the population
n_hiv_neg <- n*(1-p) ## Number of HIV-
start_year <- 2020 ## Start year of the simulation
end_year <- 2040 ## End year of the simulation
tstep <- 0.2 ## Length of a timestep (years)
nsteps <- (end_year - start_year)/tstep ## Number of time steps in the simulation
campaign_years <- seq(2020, 2040, by = 5) ## Calendar years when the campaign is implemented
campaign_tsteps <- as.integer((campaign_years - start_year)*(1/tstep) + 1) ## This just calculates the specific time steps when the campaign occurs
campaign_coverage <- 0.75 ## Proportion of target population tested by the campaign
cost_per_test <- 8.85 ## Cost ($) per person tested
n_target_pop <- n_hiv_neg + n_hiv_undiag ## Size of the target population
prop_tested_one_year <- (1-0.41) ## Proportion of incident infections that discover status within year. Based on AHRI data that 41% refuse testing within the year (you can probably find a more accurate estimate if you look harder)
## Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6340103/

## Vector of first 90 estimates to store for plotting
first_90 <- rep(NA, nsteps)
times <- start_year + tstep*((1:nsteps) - 1)

## Simulation
for(tt in 1:nsteps) {
  
  cal_year <- floor(start_year + tstep*(tt-1)) ## Calculate the current calendar year
  first_90[tt] <- prop_diagnosed ## Store this number for plotting at the end
  
  ## Print out current state of the model
  print(paste("At time step", tt, "in year", cal_year))
  print(paste("HIV prevalence =", round(p, 3)))
  print(paste("Number PLHIV diagnosed =", round(n_hiv_diag, 3)))
  print(paste("Number PLHIV undiagnosed =", round(n_hiv_undiag, 3)))
  print(paste("Proportion PLHIV diagnosed =", round(prop_diagnosed, 3)))
  
  ## Implement testing campaign (assuming equal test uptake between HIV- and undiagnosed PLHIV)
  if(tt %in% campaign_tsteps) {
    cat("\n")
    print(paste("Implementing testing campaign with", paste0((100*campaign_coverage), "%"), "coverage"))
    num_tested <- n_target_pop*campaign_coverage
    yield <- n_hiv_undiag / n_target_pop ## Percent of those tested who are undiagnosed HIV+
    print(paste("Testing", num_tested, "people with ", paste0(round(100*yield, 2), "%"), "yield at a total cost of", paste0("$", round(cost_per_test*num_tested, 2))))
    ## Update number PLHIV diagnosed based on test campaign
    n_hiv_diag <- n_hiv_diag + campaign_coverage*n_hiv_undiag
    n_hiv_undiag <- n_hiv_undiag - campaign_coverage*n_hiv_undiag
    print(paste("Campaign increased first 90 from", round(prop_diagnosed, 3), "to", round(n_hiv_diag/n_hiv, 3)))
    prop_diagnosed <- n_hiv_diag/(n_hiv_diag + n_hiv_undiag) ## Update proportion diagnosed
  }
  
  ## Simulate HIV incidence
  num_inc_infections <- inc*n_hiv_neg*tstep ## In your case, you would just calculate this from the transmission model
  print(paste(round(num_inc_infections, 3), "new infections in this time step"))
  n_hiv_undiag <- n_hiv_undiag + num_inc_infections*(1-prop_tested_one_year)
  n_hiv_diag <- n_hiv_diag + num_inc_infections*prop_tested_one_year
  print(paste("HIV incidence decreased first 90 from", round(prop_diagnosed, 3), "to", round(n_hiv_diag/(n_hiv_diag + n_hiv_undiag), 3)))
  prop_diagnosed <- n_hiv_diag/(n_hiv_diag + n_hiv_undiag) ## Update proportion diagnosed
  
  ## Recalculate population statistics (disregaring mortality, aging, etc)
  n_hiv <- n_hiv_diag + n_hiv_undiag ## Total number of PLHIV in the population
  n_hiv_neg <- n - n_hiv
  p <- n_hiv/n
  
  cat("\n\n\n")
  
}

pdf(file = "first90_time.pdf", height = 4, width = 6)
plot(x = times, y = first_90, type = "n", main = "UNAIDS First 90 over time",
     xlab = "Year", ylab = "Proportion PLHIV who are diagnosed",
     ylim = c(0, 1))
lines(x = times, y = first_90, col = "blue", lwd = 2, lty = 1)
abline(v = campaign_years, lty = "dotted")
dev.off()
