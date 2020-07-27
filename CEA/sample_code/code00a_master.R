# ***********************************************************
## GENERIC CEA CODE -----------------------------------------
# ***********************************************************

setwd("~/GD_AYA/Research Resources/CEA modeling generally/Teaching/ForMitaRianMarkAliya")

if (!dir.exists(file.path("outputs"))){
  dir.create(file.path("outputs"))
}

output_data = paste0("./outputs/output_data", "_", format(Sys.Date(), "%d%b%Y"))
output_fig_tbls = paste0("./outputs/output_fig_tbls", "_", format(Sys.Date(), "%d%b%Y"))

if (!dir.exists(file.path(output_data))){
  dir.create(file.path(output_data))
}

if (!dir.exists(file.path(output_fig_tbls))){
  dir.create(file.path(output_fig_tbls))
}

myseed = 19871219 # the seed of the analysis; this is so that every times 
# http://www.datasciencemadesimple.com/generate-sample-set-seed-function-r/
# you can put whatever seed you want. in my case, that's my birthday.
# If you want distinct threads in an MCMC, 
# this is something to be very careful with. 
# Set different seeds for each thread.
# we won't do an MCMC here. 

lbl=list()
lbl$int_names = c("No medication", "Medication I", "Medication II")

keep = c("keep", "myseed", "output_data", "output_fig_tbls", "lbl")
# It is good practice to clean up your environment periodically 
# in your workflow, BUT, there are things I want to be able to keep there.
# That's why "keep" helps us do.

# ***********************************************************
### Section 1: Set up of parameters for analysis without uncertainty ---
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 
# might be unnecessary for this section (as there is no randomness in this portion of the analysis)
# but for consistency, I set the seed here.

# B) Load output from previous sections of the analysis

# C) Run code
source("code01_pars_no_unc.R")

# D) Save output
save(fixedpars, file=paste0(output_data, "/code01_data_no_unc.Rdata"))

# E) Summary tables (NA this time)
# TODO: maybe table or tree of the fixed pars?

# ***********************************************************
### Section 2: Do CEA without uncertainty -------------------
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis
load(paste0(output_data, "/code01_data_no_unc.Rdata"))
attach(fixedpars) 
  # this makes it so that you can call the variables in that list without writing "$flexpars" over and over again 

# C) Run code
source("code02_icer.R")

# D) Save output
save(prob, ylls, ylds, dalys, cost_meds, cost_treat, costs, sum.df, ICER.df, 
     file=paste0(output_data, "/code02_cea_no_unc.Rdata"))

# E) Summary tables
# TODO: future: make a pdf of the tables with rmarkdown
# Here is a sample of the code you will need. You can see the output table in the figure window
sum.df2 = sum.df
for(i in 2:dim(sum.df)[2]){sum.df2[,i] = format(sum.df2[,i], nsmall=0, scientific= F, big.mark=",")}
tab_nounc1 = kable(sum.df2)
kable_styling(tab_nounc1, font_size = 10)
tab_nounc2 = kable(ICER.df)
kable_styling(tab_nounc2, font_size = 10)
# results of status: ND=not dominated, ED=Extendedly (Weakly) Dominated, D=Dominated.

# It is possible to make a table for a word file. Kable is really great for PDF or HTML tables, but it is less flexible for word.
# a) make it into a word file, let it be ugly (can't specify column widths), and you later edit by hand
# b) make it into an html and then copy/paste the the HTML into word. It should keep nice formating.
# c) make it into a pdf. Open the PDF with Adobe Professional or Adobe exporter and export to word. It should keep the nice formatting.
# d) To be really tailored, I recomment using flextable instead of kable:
# https://davidgohel.github.io/flextable/articles/overview.html

# ***********************************************************
### Section 3: Set up analysis with uncertainty -------------
# (net benefits framework) 
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis (NA for this section)

# C) Run code
iter= 10000
keep = c(keep, iter)
source("code03_data_unc.R")

# D) Save output
save(iter, par, file=paste0(output_data, "/code03_data_unc.Rdata"))

# E) Summary tables
# TODO: future: make a pdf of the tables with rmarkdown
tab_par_summary = kable(unlist(par_summary))
kable_styling(tab_par_summary, font_size = 10)

# ***********************************************************
### Section 4: Simulate outcomes with uncertainty -----------
# ***********************************************************
# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis
load(paste0(output_data, "/code03_data_unc.Rdata"))

# C) Run code
source("code04_sim_unc.R")

# D) Save output
save(prob, ylls, ylds, dalys, cost_meds, cost_treat, costs, sum.df,  
      file=paste0(output_data, "/code04_sim_unc.Rdata"))

# E) Summary tables
# this requires the rmarkdown package
# http://haozhu233.github.io/kableExtra/best_practice_for_newline_in_latex_table.pdf
# sum.df %>%
#   mutate_all(linebreak) %>%
#     kable("latex", booktabs = T, escape = F, align = "c")

tab1 = kable(sum.df[,1:4], escape=T)
kable_styling(tab1, font_size = 10)
tab2 = kable(sum.df[,c(1,5:7)], escape=T)
kable_styling(tab2, font_size = 10)
tab2 = kable(sum.df[,c(1,8:10)], escape=T)
kable_styling(tab2, font_size = 10)

# code for table in latex: 
# tab1 = kable(sum.df[,1:4], format = 'latex', booktabs = F, escape=T)
# kable_styling(tab1, font_size = 7)
# tab2 = kable(sum.df[,c(1,5:7)], format = 'latex', booktabs = F, escape=T)
# kable_styling(tab2, font_size = 7)
# tab2 = kable(sum.df[,c(1,8:10)], format = 'latex', booktabs = F, escape=T)
# kable_styling(tab2, font_size = 7)

# ***********************************************************
### Section 5: CEAs with uncertainty ------------------------
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis
load(paste0(output_data, "/code03_data_unc.Rdata"))
load(paste0(output_data, "/code04_sim_unc.Rdata"))

# C) Run code
source("code05_cea_unc.R")

# D) Save output
save(thresholds, nmb, ceac, ceac_pre, ceac_df, nmb_means, ceaf_prob, ceac_pie, 
     file=paste0(output_data, "/code05_cea_unc.Rdata"))

# E) Summary tables
# TODO: include EIC-style table.
# ceac_pie gives ceaC results at 0, 250, 500, 1000

# ***********************************************************
### Section 6: VOI/EVPPI ------------------------------------
# only for DISCOUNTED costs and effects
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis
load(paste0(output_data, "/code03_data_unc.Rdata"))
load(paste0(output_data, "/code04_sim_unc.Rdata"))
load(paste0(output_data, "/code05_cea_unc.Rdata"))

# C) Run code
source("code06_VOI.R")

# D) Save output
save(partial.evpi, partial.evpi_df,
     file=paste0(output_data, "/code06_VOI.Rdata"))

# E) Summary tables (NA this time)

# ***********************************************************
### Section 7: Threshold analysis ---------------------------
# ***********************************************************

# A) Set up for this section of the analysis
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
source('code00b_config.R')
set.seed(myseed) 

# B) Load output from previous sections of the analysis
load(paste0(output_data, "/code03_data_unc.Rdata"))
load(paste0(output_data, "/code04_sim_unc.Rdata"))
load(paste0(output_data, "/code05_cea_unc.Rdata"))

# C) Run code
source("code07_threshold_analysis.R")

# D) Save output (nothing in this case)

# E) Summary tables (NA this time)
