segment_plot = function(dat,title_list,filename, horiz_list){
  
  y_seq=1:dim(dat)[2]
  x_seq = seq(horiz_list$xlim[1], horiz_list$xlim[2], horiz_list$xlimsub)
  
  dat_summary = cbind(apply(dat, 2, "median"),
                      apply(dat, 2, "quantile", 0.025, na.rm=T),
                      apply(dat, 2, "quantile", 0.975, na.rm=T))
  
  bottomlabel = x_seq
  leftlabel = rev(dimnames(dat)[[2]])
  rightlabel = apply(dat_summary, 1, "ci_string_dec", ifelse(horiz_list$xlim[2]>1, 2, 4))
  rightlabel[grepl("NA", rightlabel)] = ""
  
  jpeg(filename, res=300, width=6, height=1+dim(dat)[2]/8, units="in")
  par(mar=c(2.5,6.1,2.1,6.0), oma = rep(0.5,4)) # specify margins
  plot(rev(dat_summary[,1]),y_seq, pch=20,
       ylab=title_list$ylab,xlab="",yaxt="n",xaxt="n",xlim=horiz_list$xlim,
       main=title_list$main, cex.main=0.75)
  title(xlab=title_list$xlab, line=1.1, cex.lab=0.5)
  segments(rev(dat_summary[,2]),y_seq,rev(dat_summary[,3]),y_seq)
  axis(side=2,at=y_seq,labels=leftlabel,las=1,cex.axis=0.5)
  axis(side=1,at=x_seq,labels=rep("", length(bottomlabel)),las=1,cex.axis=0.5)
  axis(side=1,at=x_seq,labels=bottomlabel ,lwd = 0, las=1,cex.axis=0.5, line=-0.75)
  axis(side=4,at=y_seq,labels=rev(rightlabel), lwd=0, las=1,cex.axis=0.5, line=-0.5)
  mtext(text="Median (95% CrI)", at=length(y_seq)+1, side=4, las=1, cex=0.5, line=0.5)
  dev.off()
  
}
