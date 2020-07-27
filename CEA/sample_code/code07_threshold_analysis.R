## Threshold analysis - probability of complication 2 ------------------------

# go back to hospital cost...
# needs: nmb and par$cost_hosp
# ceac_pre, ceaf_which could be useful
# discretize cost_hosp <400, every 25 until 700, and 700<

# put in list Ndf data.frames for each category, crunch probs and then abind

cutpts = c(seq(plyr::round_any(quantile(par$pr_comp[,2]*100, 0.05), 1),
               plyr::round_any(quantile(par$pr_comp[,2]*100, 0.95), 1), by=.5),
           plyr::round_any(quantile(par$pr_comp[,2]*100, 0.95), 1)+0.5)

test = Hmisc::cut2(par$pr_comp[,2]*100, cutpts, minmax=T)

# split up the dataset (CEAC here, CEAF below)
ta_list=lapply(1:length(unique(test)), function(j){ceac_pre[test==unique(test)[j],]})

# Now calculate summary (code is more elaborate)
# Get the ceac values for each cost_hosp category and abind them
ceac_each = function(x,y){apply(ta_list[[y]]==x, 2, "sum")/dim(ta_list[[y]])[1]}
# Get the values for ALL cost_hosp category and abind it
ceac_all = function(j){abind(lapply(1:3, function(i){ceac_each(i,j)}), along=0)}
# Finale
ceac_par = abind(lapply(1:length(unique(test)), ceac_all), along=0)

# IN one ugly line:
# ceac_par = abind(lapply(1:length(unique(test)), function(j){abind(lapply(1:3, function(x){apply(ta_list[[j]]==x, 2, "sum")/dim(ta_list[[j]])[1]}), along=0)}), along=0)

names(dim(ceac_par)) = c("par_cat", "interventions", "wtp")

dimnames(ceac_par) = list(par_cat = unique(test), 
                          interventions = lbl$int_names, 
                          wtp = thresholds)
ceac_par_df = reshape2::melt(ceac_par)

# CEAF 
taf_list = lapply(1:length(unique(test)), function(j){nmb[test==unique(test)[j],,]})

# get the mean nmb for each intervention and each wtp, stratified for par_cat
taf_pre = abind(lapply(1:length(unique(test)), function(j){apply(taf_list[[j]], 2:3, "mean")}), along=0)
taf = apply(taf_pre, c(1,3), "which.max")

dimnames(taf) = list(par_cat = unique(test), wtp = thresholds)
taf_df = reshape2::melt(taf, value.name="interventions")
taf_df$interventions = factor(taf_df$interventions)
levels(taf_df$interventions)=lbl$int_names

ceaf_tf = right_join(ceac_par_df, taf_df)

# heat map
# ggplot(ceaf_tf, aes(x=wtp, y=factor(par_cat, level = levels(unique(test))))) + 
#   themebar + theme(legend.position = "bottom") + 
#   geom_tile(aes(fill = interventions, alpha=value), color="white") + 
#   scale_alpha(range=c(0,1)) + 
#   scale_fill_manual(values=Set0) +  
#   ylab("Parameter value category") +
#   xlab("Willingness to pay \n (US$ per DALY averted)") + 
#   guides(fill = guide_legend(nrow = 1))

print(ggplot(ceaf_tf, aes(x=wtp, y=factor(par_cat, level = levels(unique(test))), colour=interventions)) + 
  themebar + theme(legend.position = "bottom") + 
  geom_point(aes(size=value)) + 
  # scale_size(range=c(0.01,1)) + 
  scale_color_manual(values=Set0) +  
  ylab("Parameter value range") +
  xlab("Willingness to pay \n (US$ per DALY averted)") + 
  labs(caption="Size of markers indicates the probability of cost-effectiveness according to ceaf") + 
  guides(fill = guide_legend(nrow = 1)))
ggsave(paste(output_fig_tbls, "/ceaf_threshold_pr_comp2.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)


# https://stackoverflow.com/questions/13016022/ggplot2-heatmaps-using-different-gradients-for-categories

## Threshold analysis: Probability of complication 3 -----------------

# go back to yll...
# needs: nmb and par$cost_hosp
# ceac_pre, ceaf_which could be useful
# discretize cost_hosp <400, every 25 until 700, and 700<

# put in list Ndf data.frames for each category, crunch probs and then abind

cutpts = c(seq(plyr::round_any(quantile(par$pr_comp[,3]*100, 0.05), 0.5),
               plyr::round_any(quantile(par$pr_comp[,3]*100, 0.95), 0.5), by=0.5),
           plyr::round_any(quantile(par$pr_comp[,3]*100, 0.95), 0.5)+0.5)

test = Hmisc::cut2(par$pr_comp[,3]*100, cutpts, minmax=T)

# split up the dataset (CEAC here, CEAF below)
ta_list=lapply(1:length(unique(test)), function(j){ceac_pre[test==unique(test)[j],]})

# Now calculate summary (code is more elaborate)
# Get the ceac values for each cost_hosp category and abind them
ceac_each = function(x,y){apply(ta_list[[y]]==x, 2, "sum")/dim(ta_list[[y]])[1]}
# Get the values for ALL cost_hosp category and abind it
ceac_all = function(j){abind(lapply(1:3, function(i){ceac_each(i,j)}), along=0)}
# Finale
ceac_par = abind(lapply(1:length(unique(test)), ceac_all), along=0)

# IN one ugly line:
# ceac_hosp = abind(lapply(1:length(unique(test)), function(j){abind(lapply(1:3, function(x){apply(ta_list[[j]]==x, 2, "sum")/dim(ta_list[[j]])[1]}), along=0)}), along=0)

names(dim(ceac_par)) = c("par_cat", "interventions", "wtp")
dimnames(ceac_par) = list(par_cat = unique(test), 
                          interventions = lbl$int_names, 
                          wtp = thresholds)
ceac_par_df = reshape2::melt(ceac_par)

# CEAF 
taf_list = lapply(1:length(unique(test)), function(j){nmb[test==unique(test)[j],,]})

# get the mean nmb for each intervention and each wtp, stratified for hosp_cost
taf_pre = abind(lapply(1:length(unique(test)), function(j){apply(taf_list[[j]], 2:3, "mean")}), along=0)
taf = apply(taf_pre, c(1,3), "which.max")

dimnames(taf) = list(par_cat = unique(test), wtp = thresholds)
taf_df = reshape2::melt(taf, value.name="interventions")
taf_df$interventions = factor(taf_df$interventions)
levels(taf_df$interventions)=lbl$int_names

ceaf_tf = right_join(ceac_par_df, taf_df)

# heat map
# ggplot(ceaf_tf, aes(x=wtp, y=factor(par_cat, level = levels(unique(test))))) + 
#   themebar + theme(legend.position = "bottom") + 
#   geom_tile(aes(fill = interventions, alpha=value), color="white") + 
#   scale_alpha(range=c(0,1)) + 
#   scale_fill_manual(values=Set0) +  
#   ylab("YLL category") +
#   xlab("Willingness to pay \n (US$ per DALY averted)") + 
#   guides(fill = guide_legend(nrow = 1))

print(ggplot(ceaf_tf, aes(x=wtp, y=factor(par_cat, level = levels(unique(test))), colour=interventions)) + 
        themebar + theme(legend.position = "bottom") + 
        geom_point(aes(size=value)) + 
        # scale_size(range=c(0.01,1)) + 
        scale_color_manual(values=Set0) +  
        ylab("Parameter value range") +
        xlab("Willingness to pay \n (US$ per DALY averted)") + 
        labs(caption="Size of markers indicates the probability of cost-effectiveness according to ceaf") + 
        guides(fill = guide_legend(nrow = 1)))
ggsave(paste(output_fig_tbls, "/ceaf_threshold_pr_comp3.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

# https://stackoverflow.com/questions/13016022/ggplot2-heatmaps-using-different-gradients-for-categories

## Threshold YLL -------------------------------------------

# Threshold analysis

# go back to yll...
# needs: nmb and par$cost_hosp
# ceac_pre, ceaf_which could be useful
# discretize cost_hosp <400, every 25 until 700, and 700<

# put in list Ndf data.frames for each category, crunch probs and then abind

cutpts = c(seq(plyr::round_any(quantile(par$yll, 0.05), 0.5),
               plyr::round_any(quantile(par$yll, 0.95), 0.5), by=0.5),
           plyr::round_any(quantile(par$yll, 0.95), 0.5)+0.5)

test = Hmisc::cut2(par$yll, cutpts, minmax=T)

# split up the dataset (CEAC here, CEAF below)
ta_list=lapply(1:length(unique(test)), function(j){ceac_pre[test==unique(test)[j],]})

# Now calculate summary (code is more elaborate)
# Get the ceac values for each cost_hosp category and abind them
ceac_each = function(x,y){apply(ta_list[[y]]==x, 2, "sum")/dim(ta_list[[y]])[1]}
# Get the values for ALL cost_hosp category and abind it
ceac_all = function(j){abind(lapply(1:3, function(i){ceac_each(i,j)}), along=0)}
# Finale
ceac_par = abind(lapply(1:length(unique(test)), ceac_all), along=0)

# IN one ugly line:
# ceac_hosp = abind(lapply(1:length(unique(test)), function(j){abind(lapply(1:3, function(x){apply(ta_list[[j]]==x, 2, "sum")/dim(ta_list[[j]])[1]}), along=0)}), along=0)

names(dim(ceac_par)) = c("yll_cat", "interventions", "wtp")
dimnames(ceac_par) = list(yll_cat = unique(test), 
                          interventions = lbl$int_names, 
                          wtp = thresholds)
ceac_par_df = reshape2::melt(ceac_par)

# CEAF 
taf_list = lapply(1:length(unique(test)), function(j){nmb[test==unique(test)[j],,]})

# get the mean nmb for each intervention and each wtp, stratified for hosp_cost
taf_pre = abind(lapply(1:length(unique(test)), function(j){apply(taf_list[[j]], 2:3, "mean")}), along=0)
taf = apply(taf_pre, c(1,3), "which.max")

dimnames(taf) = list(yll_cat = unique(test), wtp = thresholds)
taf_df = reshape2::melt(taf, value.name="interventions")
taf_df$interventions = factor(taf_df$interventions)
levels(taf_df$interventions)=lbl$int_names

ceaf_tf = right_join(ceac_par_df, taf_df)

# heat map
# ggplot(ceaf_tf, aes(x=wtp, y=factor(yll_cat, level = levels(unique(test))))) + 
#   themebar + theme(legend.position = "bottom") + 
#   geom_tile(aes(fill = interventions, alpha=value), color="white") + 
#   scale_alpha(range=c(0,1)) + 
#   scale_fill_manual(values=Set0) +  
#   ylab("YLL category") +
#   xlab("Willingness to pay \n (US$ per DALY averted)") + 
#   guides(fill = guide_legend(nrow = 1))

print(ggplot(ceaf_tf, aes(x=wtp, y=factor(yll_cat, level = levels(unique(test))), colour=interventions)) + 
        themebar + theme(legend.position = "bottom") + 
        geom_point(aes(size=value)) + 
        # scale_size(range=c(0.01,1)) + 
        scale_color_manual(values=Set0) +  
        ylab("Parameter value range") +
        xlab("Willingness to pay \n (US$ per DALY averted)") + 
        labs(caption="Size of markers indicates the probability of cost-effectiveness according to ceaf") + 
        guides(fill = guide_legend(nrow = 1)))
ggsave(paste(output_fig_tbls, "/ceaf_threshold_yll.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

# https://stackoverflow.com/questions/13016022/ggplot2-heatmaps-using-different-gradients-for-categories

