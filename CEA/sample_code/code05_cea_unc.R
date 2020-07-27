## CEA with uncertainty ---------------------------

## CE-plane graph ---------------------------------

costs_total = apply(sweep(prob*costs, 1, par$sick, "*"), c(1,3), "sum", na.rm=T)
dalys_total = apply(sweep(prob*dalys, 1, par$sick, "*"), c(1,3), "sum", na.rm=T)

dalys_ceplane = sweep(-dalys_total, 1, dalys_total[,1], "+")
costs_ceplane = sweep(costs_total, 1, costs_total[,1], "-")

# sanity check: #1 in the second dimension should be all zeros
# summary(dalys_ceplane[,1])
# summary(costs_ceplane[,1])

# Table of differences (costs, dalys, cases, deaths)?

# ggplot for CE-plane (didactic...)
# melt and join arrays then feed to ggplot

ceplane = inner_join(reshape2::melt(costs_ceplane, value.name="cost_dif"), 
                     reshape2::melt(dalys_ceplane, value.name="daly_dif"))

# NOTE: you need to put the ggplot command in a print() command so that it is made when you call this code within source()
png(paste(output_fig_tbls, "/cea_plane_traditional.png",sep=""), 
    width = 7, height = 5, units="in", res=300)
print(ggplot(data=ceplane[ceplane$iteration<300,], aes(x=daly_dif, y=cost_dif/1000, colour=interventions)) + #, colour=intervention
  geom_point(size=1.5, alpha=0.3, stroke=0.5) + # facet_wrap(~site, scales="free", ncol=1) + #
  themebar + theme(legend.position = "bottom") + 
  ylab("Difference in costs, thousands $") +
  xlab("Disability-adjusted life-years (DALYs) averted") + 
  scale_y_continuous(labels=scales::dollar) +
  scale_x_continuous(labels=scales::comma) +
  scale_colour_manual(values=Set0, guide=guide_legend(nrow=1)))
dev.off()

cetplane = inner_join(reshape2::melt(costs_total, value.name="cost"), 
                      reshape2::melt(dalys_total, value.name="daly"))

png(paste(output_fig_tbls, "/cea_plane_total.png",sep=""), 
    width = 7, height = 5, units="in", res=300)
print(ggplot(data=cetplane[cetplane$iteration<300,], aes(x=daly, y=cost/1000, colour=interventions)) + #, colour=intervention
  geom_point(size=1.5, alpha=0.3, stroke=0.5) + # facet_wrap(~site, scales="free", ncol=1) + #
  themebar + theme(legend.position = "bottom") + 
  ylab("Difference in costs, thousands $") +
  xlab("Disability-adjusted life-years (DALYs) averted") + 
  scale_y_continuous(labels=scales::dollar) +
  scale_x_continuous(labels=scales::comma) +
  scale_colour_manual(values=Set0, guide=guide_legend(nrow=1)))
dev.off()
# apply(costs_ceplane, 2, "quantile", c(0.5, 0.025, 0.975))
# apply(dalys_ceplane, 2, "quantile", c(0.5, 0.025, 0.975))

## CEAC -------------------------------------

thresholds = seq(0, 1000, 50)
ceac_array_0 = array(0, dim=c(iterations = iter, 
                              interventions = length(lbl$int_names), 
                              thresholds = length(thresholds)), 
                     dimnames = list(iterations=NULL, 
                                     interventions=lbl$int_names,
                                     thresholds = thresholds))

nmb = sweep(sweep(ceac_array_0, 1:2, dalys_ceplane, "+"), 3, thresholds, "*") - 
  sweep(ceac_array_0, 1:2, costs_ceplane, "+")

ceac_pre = apply(nmb, c(1,3), "which.max")
ceac = abind(lapply(1:3, function(x){apply(ceac_pre==x, 2, "sum")}), along=0)/iter

# label ceac
names(dimnames(ceac)) = c("interventions", "wtp")
dimnames(ceac) = list(interventions=lbl$int_names, wtp=thresholds)
# melt
ceac_df = reshape::melt(ceac)
ceac_df$interventions = factor(ceac_df$interventions, 
                               levels=levels(ceac_df$interventions)[c(3,1,2)])

print(ggplot(data=ceac_df, aes(x=wtp, y=value, colour=interventions)) +
  geom_line(alpha=1) +
  themebar + theme(legend.position = "bottom") +
  scale_colour_manual(values=Set0) +
  guides(colour = guide_legend(nrow = 1)) +
  ylab("Probability Cost-effective") + scale_y_continuous(limits=c(0,1)) +
  xlab("Willingness to pay (US$ per DALY averted)"))
ggsave(paste(output_fig_tbls, "/ceac_traditional.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

# CEAF -----------------------------------
# for each WTP (and for each scenario), which of the 
# interventions has the largest *mean* NB? Then plot.

# back to using nmb
nmb_means = apply(nmb, 2:length(dim(nmb)), "mean")

# sanity check
# matplot(t(nmb_means), type="l")

ceaf_which = apply(nmb_means, 2, "which.max")

ceaf_which_df = data.frame(interventions=ceaf_which, wtp=thresholds)
ceaf_which_df$interventions = factor(ceaf_which_df$interventions)
levels(ceaf_which_df$interventions)=lbl$int_names

ceaf_prob = right_join(ceac_df, ceaf_which_df)

# Plots CEAC, CEAF ---------------------
# TODO: Add CEAC vs CEAF to the legend.
# http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-softwarea
# TODO: Make additional graph with horizontal lines to mark what pies depict

print(ggplot(data=ceac_df, aes(x=wtp, y=value, colour=interventions)) +
        geom_line(alpha=1) + 
        themebar + theme(legend.position = "bottom") + 
        scale_colour_manual(values=Set0) +  
        guides(colour = guide_legend(nrow = 1)) + 
        geom_line(data=ceaf_prob, aes(x=wtp, y=value*1.02), color="black", linetype="dashed") + 
        ylab("Probability cost-effective") + scale_y_continuous(limits=c(0,1.05)) + 
        xlab("Willingness to pay \n (US$ per DALY averted)"))
ggsave(paste(output_fig_tbls, "/ceaf_traditional.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

## Pies ----------------------------
# this is something we like in HATMEPP, but it is not traditional

ceac_pie = ceac_df[ceac_df$wtp %in% c("0", "250", "500", "1000"),]

# for position of text
# https://stackoverflow.com/questions/24803460/r-ggplot2-add-labels-on-facet-pie-chart

ceac_pie$wtp = paste("Willingness to pay:\n$", ceac_pie$wtp, " per DALY averted", sep="")

ceac_pie$wtp = factor(ceac_pie$wtp, levels=unique(ceac_pie$wtp))

print(ggplot(data=ceac_pie, 
       aes(x=factor(1), y=value, fill=interventions, group=interventions)) +
  geom_bar(stat="identity") +  
  themebar + theme(legend.position = "bottom", 
                   axis.text.x=element_blank(), 
                   axis.text.y=element_blank(), 
                   axis.ticks.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.title.x=element_blank(),
                   strip.text = element_text(size=rel(0.75), face="bold"),
                   legend.text = element_text(size=rel(0.75), face="bold")) + 
  scale_fill_manual(values=Set0) +  
  guides(fill = guide_legend(nrow = 1)) + 
  facet_wrap(.~wtp, nrow=2) + 
  coord_polar("y", start=0))
# geom_text(aes(y = value/2, 
#     label = scales::percent(value))) + 
ggsave(paste(output_fig_tbls, "/ceac_pies.eps",sep=""), 
       width = 7, height = 5, units="in", dpi=400)

