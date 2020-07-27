# Set up of parameters ------------------------------------------

# change parameter values for each iter
par=list()

tmp0 = find_gamma_pars(sick/2, sick*2)
par$sick =  rgamma(iter, tmp0[1], scale=tmp0[2])

tmp_beta_samp1 = rbeta(iter, pr_comp[1]*100, 100-pr_comp[1]*100)
tmp_beta_samp2 = rbeta(iter, 90, 10)
tmp_beta_samp3 = rnorm(iter, 0.975, 0.025)

par$pr_comp = cbind(tmp_beta_samp1, 
                    tmp_beta_samp1*tmp_beta_samp2,
                    tmp_beta_samp1*tmp_beta_samp3)

tmp_beta_samp1 = rbeta(iter, pr_death[1]*100+1, 100-pr_death[1]*100+1)
tmp_beta_samp2 = rnorm(iter, 0.975, 0.025)
tmp_beta_samp3 = rbeta(iter, 80, 20)

par$pr_death = cbind(tmp_beta_samp1, 
                     tmp_beta_samp1*tmp_beta_samp2,
                     tmp_beta_samp1*tmp_beta_samp3)

par$pr_recover = 1-par$pr_death

tmp1 = find_gamma_pars(yld_comp/2, yld_comp*2)
par$yld_comp = rgamma(iter, tmp1[1], scale=tmp1[2])

par$yld_simple = par$yld_comp*rbeta(iter, 50, 50)

par$yll = rgamma(iter, yll^2/((yll/4)^2), scale=((yll/4)^2)/yll) 

# if you want to parameterize a gamma distribution with cis, 
# use find_gamma_pars or find_gamma_flex from code00b_functions
tmp2 = find_gamma_pars(cost_hosp/3, cost_hosp*3) 
par$cost_hosp = rgamma(iter, tmp2[1], scale=tmp2[2])

par$cost_rx = cost_rx

# Sanity check -----------------
# check what the statistical summaries are of all the parameters are 
tmp_fcn = function(par_list){
  if(is.vector(par_list)){tmpobj = data.frame(par_list)}else{tmpobj = par_list}
  tmpchar = paste("Mean: ", format(apply(tmpobj, 2, "mean"), digits=2, nsmall = 2, scientific=F), ", Median (95p CI): ",
                  apply(apply(tmpobj, 2, "quantile", c(0.5, 0.025, 0.975)), 2, "ci_string_dec", 2), sep="")
  return(tmpchar)}

par_summary = lapply(par, "tmp_fcn")


