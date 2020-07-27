dist_hist_panel = function(list_dist, desc_df, store_sc){
  # list_dist is a list with the vectors of the distribution
  # desc_df is a data.frame with the descriptions of the paraters in list_dist
  # store_rc is the destination of the histograms
jpeg(paste(store_sc, unique(desc_df$par_grp[1]), ".jpg", sep=""), res=300, width=10, height=10, units="in")
par(mfrow=c(floor(sqrt(dim(desc_df)[1])), ceiling(dim(desc_df)[1]/floor(sqrt(dim(desc_df)[1])))), mar=c(2,2,5,2))
for (i in 1:dim(desc_df)[1]){
  tmp_vector = list_dist[[desc_df$var_label[i]]]
  
  par_summ1 = paste("95% CI: ", format(median(tmp_vector), digits=2, nsmall=2), " (",
                    format(quantile(tmp_vector, 0.025), digits=2, nsmall=2), ", ", 
                    format(quantile(tmp_vector, 0.975), digits=2, nsmall=2), ")", sep="")
  par_summ2 = paste("Min: ",
                    format(min(tmp_vector), digits=2, nsmall=2), " Max: ", 
                    format(max(tmp_vector), digits=2, nsmall=2), sep="")
  dist_note = paste(desc_df$rcode_dist[i], "(", desc_df$rcode_par1[i], ",", 
                    desc_df$rcode_par2[i], ")", sep="")
  
  hist(tmp_vector, xlim=c(quantile(tmp_vector, 0.025), quantile(tmp_vector, 0.975)), breaks=100, 
       freq=F, col="darkblue", xlab="", main="", ylab="", cex.axis=0.7) 
  
  mtext(wrapper(desc_df$var_desc[i], 60), 3, line = 3.25, cex=0.5, font=2)
  mtext(desc_df$var_label[i], 3, line = 2.25, cex=0.5, font=3)
  mtext(dist_note, 3, line = 1.5, cex=0.5, font=1)
  mtext(par_summ1, 3, line = 0.75, cex=0.5, font=1)
  mtext(par_summ2, 3, line = 0, cex=0.5, font=1)
}
dev.off()
}