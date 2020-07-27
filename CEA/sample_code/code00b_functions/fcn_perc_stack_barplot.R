perc_stack_barplot = function(dat,title_list,filename){
  
  # dat: 2-dimensional data to plot
  # title_list: list; must contain main title, ylab, xlab
  # filename: file name (including path)

jpeg(filename, res=300, width=6, height=4, units="in")
par(xpd=TRUE, mar=c(4.5,3.1,1.1,1.1))
barplot(dat, col=2:(dim(dat)[1]+1), border=NA, cex.axis=0.5, 
        main=title_list$main, cex.main=0.5, cex.lab=0.5) 
title(ylab=title_list$ylab, line=1.5, cex.lab=0.5)
title(xlab=title_list$xlab, line=0, cex.lab=0.5)
legend(0, -.1*max(colSums(dat)), dimnames(dat)[[1]], col=2:(dim(dat)[1]+1), pch=20, ncol=4, cex=.5)
dev.off()

}