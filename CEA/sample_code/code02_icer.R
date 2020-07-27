# make arrays: prob, ylls, ylds, costs_meds, costs_treatment -------------

# where is the intervention dimension?
prob = array(NA, dim=c(outcome=3, interventions=length(lbl$int_names)), 
             dimnames=list(outcome=c("Recover_Comp", "Death", "Recover_Simp"),
                           interventions=lbl$int_names))
ylls=ylds=cost_meds=cost_treat=prob

# Calculate the values that go in the array ------------------------------
for(i in 1:length(lbl$int_names)){prob[,i]=c(pr_comp[i]*pr_recover[i], pr_comp[i]*(1-pr_recover[i]), 1-pr_comp[i])}
for(i in 1:length(lbl$int_names)){ylds[,i]=c(yld_comp, yld_comp, yld_simple)}
for(i in 1:length(lbl$int_names)){ylls[,i]=c(0, yll, 0)}
for(i in 1:length(lbl$int_names)){cost_treat[,i]=c(cost_hosp, cost_hosp, 0)}
for(i in 1:length(lbl$int_names)){cost_meds[,i]=rep(cost_rx[i], 3)}

# Two more arrays that are sums of other arrays
dalys = ylds + ylls
costs = cost_treat+cost_meds

# Put it all in a nice table that can be displayed -----------------------
sum.df = data.frame(Interventions = lbl$int_names)
sum.df$`Patients` = sick
sum.df$`Hospitalised` = sick*(1-prob[3,])
sum.df$`Deaths` = sick*prob[2,]
sum.df$`Treat $` = apply(sick*prob*cost_treat, 2, "sum", na.rm=T)
sum.df$`Med $` = apply(sick*prob*cost_meds, 2, "sum", na.rm=T)
sum.df$`Total $` = apply(sick*prob*costs, 2, "sum", na.rm=T)
sum.df$`YLLs` = apply(sick*prob*ylls, 2, "sum", na.rm=T)
sum.df$`YLDs` = apply(sick*prob*ylds, 2, "sum", na.rm=T)
sum.df$`DALYs` = apply(sick*prob*dalys, 2, "sum", na.rm=T)

# Calculate ICERs --------------------------------------------------------
ICER.df = data.frame(Strategy = c("Comparator", "Med I", "Med II"),
                     Cost = c(0, sum.df$`Total $`[2:length(lbl$int_names)]-sum.df$`Total $`[1]), 
                     Effect = c(0, sum.df$DALYs[1]-sum.df$DALYs[2:length(lbl$int_names)]))

tmp_icers1 = compute_icers(ICER.df)
tmp_icers2 = calculate_icers(ICER.df)
# results of status: ND=not dominated, ED=Extendedly (Weakly) Dominated, D=Dominated.

ICER.df$ACER = ICER.df$Cost/ICER.df$Effect
ICER.df = full_join(ICER.df, tmp_icers2[,c(1,4,5)])
