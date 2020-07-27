# Configurations for analysis ---------------------------

# ********************************************************
# Necessary packages -------------------------------------
# ********************************************************

# command below: what repo are we using? 
# getOption("repos")
# Repos occasionally have glitches. One quick attempt at 
# trouble shooting is changing the repo.

library(RColorBrewer)
library(scales)
library(grid)
# library(tidyverse)
# installed as part of tidyverse: 
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)

library(rmarkdown)
library(knitr) 
library(kableExtra) # kable comes within the knitr library
library(stringr)
library(flextable)

library(readxl)
library(openxlsx)
library(R.matlab)

library(reshape2) # DEPRECATED
library(data.table)

library(MASS) # for fitdistr()
library(mvtnorm)
library(meta)
library(metafor)
library(Hmisc)
library(mgcv)

library(abind)
library(fmsb) # for radarchart (used for sanity checks only)

library(citr) 

# these should be in the packrat library of the project
# https://rstudio.github.io/packrat/walkthrough.html

# ********************************************************
# USEFUL FUNCTIONS --------------------------------------
# ********************************************************

# Logistic transformation
logit = function(x){
  x[x<0.001] = 0.001
  x[x>0.999] = 0.999
  log(x/(1-x))
}
logistic = function(x){exp(x)/(1+exp(x))}
logit(c(0, 0.5, 1)) # sanity check

# Round up to nearest specified number (useful to specify plot axes)
roundUp = function(x,to=10){to*(x%/%to + as.logical(x%%to))}

# this is useful to melt with tidyr
namedim1=function(x){
  dimnames(x)[[1]] = c(iterations=as.character(1:dim(x)[1]))
  return(x)}

fls = list.files("./code00b_functions/", pattern="^[fcn]")
for (i in 1:length(fls)){source(paste0("./code00b_functions/", fls[i]))}

# The cea-calculator is now a function, and 
# the file has the citation/link to the original.
source("./code00b_functions/icer_calculator-master/icer_calculator.R")

# ********************************************************
# GGPLOT settings ---------------------------------------
# ********************************************************

# Graphic specifications for ggplot:
themebar = theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
                 axis.title.x = element_text(size = 10, angle = 0, face="bold"),
                 axis.text.y = element_text(face="bold", color="black", size=10, angle=0), 
                 axis.title.y = element_text(size = 10, angle = 90, face="bold"),
                 plot.title = element_text(size = 16, angle = 0, face="bold"),
                 plot.subtitle = element_text(size = 14, angle = 0, face="bold"),
                 panel.border = element_rect(linetype = "solid", colour = "black", fill=NA),
                 legend.text = element_text(size = 12, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm')),
                 legend.position = "bottom",
                 legend.box = "vertical",
                 legend.background = element_rect(fill=NA, size=0.25, linetype="solid", colour ="black"),
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 legend.spacing.x=unit(0,"cm"), 
                 panel.grid.major = element_line(colour="gray", linetype = "dotted"),
                 panel.background = element_rect(fill = NA),
                 strip.background = element_rect(fill = NA),
                 strip.text = element_text(size=12, face="bold")) 

# Colors for the interventions -------------------------
library("RColorBrewer")

# display.brewer.pal(n = 9, name = 'Set1') # Only if I want to visualize them
Set0 = brewer.pal(n = 9, name = "Set1") # Include do-nothing option, which will be delineated in gray.
Set0 = Set0[c(9, 1:8)]

Set1 = brewer.pal(n = 9, name = "Set1")

