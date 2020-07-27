# Value of information ---------------------------------

# make par into dataframe ------------------------------
# do not keep par$pr_death and par$pr_recover
uncpar=which(!(names(par) %in% c("pr_recover", "cost_rx")))

par.df = data.frame(par[[uncpar[1]]])
if(!is.null(dim(par[[uncpar[1]]])[1])){
  colnames(par.df)[colnames(par.df)==""] = paste(names(par)[uncpar[1]], 1:dim(par[[uncpar[j]]])[2], sep="")
} else {
  colnames(par.df) = names(par)[uncpar[1]]                 
}

for(j in 2:length(uncpar)){
  tmp_size = dim(par.df)[2]
  par.df = cbind(par.df, par[[uncpar[j]]])
  if(!is.null(dim(par[[uncpar[j]]])[1])){
    colnames(par.df)[(tmp_size+1):(dim(par.df)[2])] = paste(names(par)[uncpar[j]], 1:dim(par[[uncpar[j]]])[2], sep="")
  } else {
    colnames(par.df)[tmp_size+1] = names(par)[uncpar[j]]                 
  }
}

# g.hat_new = array(0, c(dim(nmb)[1:5], length(lambdas), dim(unc_pars_keep)[2]))
g.hat_new = array(0, c(iter, dim(nmb)[2], length(thresholds), dim(par.df)[2]))

b = Sys.time()
# loop for parameters
library(gamm4)

for(l in 2:dim(nmb)[2]){ # for each of the interventions (that are not the referent)
  for(m in 1:length(uncpar)){ # for each parameter
    for(n in 1:dim(nmb)[3]){ # for each wtp
      g.hat_new[,l,n,m] = gam(nmb[,l,n] ~ par.df[,m])$fitted
    }
  }
}

a = Sys.time()
# a-b 

perfect.info = apply(apply(g.hat_new,c(1,3,4), "max"), c(2,3), "mean")
baseline = apply(apply(g.hat_new, 2:4,"mean"), c(2,3), "max")

partial.evpi = perfect.info - baseline ## estimate EVPI 

## Graph -----------------------------------------------

# with base plot:
# par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE)
# 
# matplot(thresholds, partial.evpi, type="l")
# legend(1200, 30000, legend=colnames(par.df), col=1:11,lty=1:11,
#        title="Ten Parameters with Highest EVPPI", cex=0.5)

# with ggplot:
dimnames(partial.evpi) = list(wtp = thresholds, pars=colnames(par.df))
partial.evpi_df = reshape2::melt(partial.evpi, value.name="partial_evpi")

print(ggplot(data=partial.evpi_df, aes(x=wtp, y=partial_evpi, colour=pars)) +
  geom_line() + # aes(linetype=pars)
  themebar + theme(legend.position = "right") + 
  scale_colour_manual(values=colorRamps::primary.colors(12, steps=4)) +  
  guides(colour = guide_legend(ncol = 1)) + 
  ylab("Partial EVPI") + # scale_y_continuous(limits=c(0,1.05)) + 
  xlab("Willingness to pay \n (US$ per DALY averted) \n Threshold"))
ggsave(paste(output_fig_tbls, "/ceac_pies.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

