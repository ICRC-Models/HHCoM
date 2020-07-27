
# Simulation arrays ------------------------------------
prob = array(NA, dim=c(iterations = iter, outcome=3, interventions=length(lbl$int_names)),
             dimnames=list(iterations=NULL, outcome=c("Recover_Comp", "Death", "Recover_Simp"),
                           interventions=lbl$int_names))
ylls=ylds=cost_meds=cost_treat=prob

for (j in 1:iter){
  
  unc.sick = par$sick[j]
  unc.pr_comp = par$pr_comp[j, ] # one for each intervention
  unc.pr_recover = par$pr_recover[j, ] # one for each intervention
  unc.yld_comp = par$yld_comp[j]
  unc.yld_simple = par$yld_simple[j]
  unc.yll = par$yll[j]
  unc.cost_hosp = par$cost_hosp[j]
  # cost_rx # one for each intervention - no uncertainty.
  
  for(i in 1:length(lbl$int_names)){prob[j,,i]=c(unc.pr_comp[i]*unc.pr_recover[i], unc.pr_comp[i]*(1-unc.pr_recover[i]), 1-unc.pr_comp[i])}
  for(i in 1:length(lbl$int_names)){ylds[j,,i]=c(unc.yld_comp, unc.yld_comp, unc.yld_simple)}
  for(i in 1:length(lbl$int_names)){ylls[j,,i]=c(0, unc.yll, 0)}
  for(i in 1:length(lbl$int_names)){cost_treat[j,,i]=c(unc.cost_hosp, unc.cost_hosp, 0)}
  for(i in 1:length(lbl$int_names)){cost_meds[j,,i]=rep(cost_rx[i], 3)}
  
}

dalys = ylds + ylls
costs = cost_treat+cost_meds

# summaries -----------------------------
# this requires the nice comma print functions from the code00b_functions
sum.df = data.frame(Interventions = lbl$int_names)
sum.df$`Sick` = apply(apply(replicate(length(lbl$int_names), par$sick), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`Hospitalised` = apply(apply(sweep(1-prob[,3,], 1, par$sick, "*"), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`Deaths` = apply(apply(sweep(prob[,2,], 1, par$sick, "*"), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`Treat $` = apply(apply(apply(sweep(prob*cost_treat, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`Medication $` = apply(apply(apply(sweep(prob*cost_meds, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`Total $` = apply(apply(apply(sweep(prob*costs, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`YLLs` = apply(apply(apply(sweep(prob*ylls, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`YLDs` = apply(apply(apply(sweep(prob*ylds, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
sum.df$`DALYs` = apply(apply(apply(sweep(prob*dalys, 1, par$sick, "*"), c(1,3), "sum", na.rm=T), 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_comma")
