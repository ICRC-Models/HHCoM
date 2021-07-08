### Script to process output on one-time screening optimal age ###

# This script takes output from the simulation model and processes it to produce plots and statistics for tables and text. 
# The script expects four output files from each scenario:
# 1. Cervical cancer cases
# 2. Cervical cancer incidence rates (raw)
# 3. Cervical cancer incidence rates (age-standardized)
# 4. HIV prevalence in females Cervical cancer 
# For each, there should be separate tabs/sheets for each of the following groups: the Total population, HIV-negative women, HIV-positive untreated women, HIV-positive women on ART, all HIv-positive women __IN THAT ORDER__. 
# Within each tab column 1 should show the year of simulation from 2020 to 2070. Column 2 should show the outcome for the relevant year, group, and scenario

## The filepaths and names correspond to these scenarios:
# S0 = Cytology at age 35 with 48% coverage throughout the timeframe. 9v HPV vaccination at 57%
# S0b = Cytology at ages 35 and 45 with 48% coverage throughout the timeframe. 9v HPV vaccination at 57%
# S1 = HPV screening at ages 35 and 45 scaled up to 90% by 2045 throughout the timeframe. 9v HPV vaccination at 90%
# S2 = S2 but with screening for WLHIV every 3 years ages 25-49
# S3 = S2 but with catch-up vaccination for WLHIV ages 15-24 at 50%

##########
# SETUP # ----
##########

## Load packages (if the packages are not yet installed, type: install.packages("packagename"))
library("ggplot2")
library("tidyverse")
library("knitr")
library("kableExtra")
library("readxl")
library("reshape2")
library("colorspace")
library("viridis")


##########
# SPECIFY FILEPATHS # ----
##########

## Set directory to folder where output is stored <<--
# setwd("C:/Users/dpwhite/Dropbox/HPV and cervical cancer modeling/Papers/PAF from HIV/Model output")
# For use on my Mac
setwd("/Users/darcywhite/Dropbox/HPV and cervical cancer modeling/Papers/PAF from HIV/Model output")

## Set filepaths <<--
# List file names for files containing incidence output for each scenario
# These files are expected to have one tab for each population subgroup (All, HIV-, untreated HIV+, treated HIV+, all HIV+)
# and for each tab to list the years in column 1 and the outcome in column 2
S0.stART.inc.path <- "PAF_crudeAnnualCC_S0.xlsx"
S0.stART.ir.path <- "PAF_crudeICC_S0.xlsx"
S0.stART.sir.path <- "PAF_ASICC_S0.xlsx"
S0.stART.hiv.path <- "PAF_crudeHivFaged15plus_S0.xlsx"

S1.stART.inc.path <- "PAF_crudeAnnualCC_S1.xlsx"
S1.stART.ir.path <- "PAF_crudeICC_S1.xlsx"
S1.stART.sir.path <- "PAF_ASICC_S1.xlsx"
S1.stART.hiv.path <- "PAF_crudeHivFaged15plus_S1.xlsx"

S2.stART.inc.path <- "PAF_crudeAnnualCC_S2.xlsx"
S2.stART.ir.path <- "PAF_crudeICC_S2.xlsx"
S2.stART.sir.path <- "PAF_ASICC_S2.xlsx"
S2.stART.hiv.path <- "PAF_crudeHivFaged15plus_S2.xlsx"

S3.stART.inc.path <- "PAF_crudeAnnualCC_S3.xlsx"
S3.stART.ir.path <- "PAF_crudeICC_S3.xlsx"
S3.stART.sir.path <- "PAF_ASICC_S3.xlsx"
S3.stART.hiv.path <- "PAF_crudeHivFaged15plus_S3.xlsx"

S4.stART.inc.path <- "PAF_crudeAnnualCC_S4.xlsx"
S4.stART.ir.path <- "PAF_crudeICC_S4.xlsx"
S4.stART.sir.path <- "PAF_ASICC_S4.xlsx"
S4.stART.hiv.path <- "PAF_crudeHivFaged15plus_S4.xlsx"

## Specify whether the excel sheets contain a header row with variable names (TRUE/FALSE)  <<--
header <- FALSE

##########
# READ AND PROCESS MODEL OUTPUT FILES # ----
##########

## Load model output files
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, col_names = header))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

S0.stART.inc <- read_excel_allsheets(S0.stART.inc.path)
S0.stART.ir <- read_excel_allsheets(S0.stART.ir.path)
S0.stART.sir <- read_excel_allsheets(S0.stART.sir.path)
S0.stART.hiv <- read_excel_allsheets(S0.stART.hiv.path)

S1.stART.inc <- read_excel_allsheets(S1.stART.inc.path)
S1.stART.ir <- read_excel_allsheets(S1.stART.ir.path)
S1.stART.sir <- read_excel_allsheets(S1.stART.sir.path)
S1.stART.hiv <- read_excel_allsheets(S1.stART.hiv.path)

S2.stART.inc <- read_excel_allsheets(S2.stART.inc.path)
S2.stART.ir <- read_excel_allsheets(S2.stART.ir.path)
S2.stART.sir <- read_excel_allsheets(S2.stART.sir.path)
S2.stART.hiv <- read_excel_allsheets(S2.stART.hiv.path)

S3.stART.inc <- read_excel_allsheets(S3.stART.inc.path)
S3.stART.ir <- read_excel_allsheets(S3.stART.ir.path)
S3.stART.sir <- read_excel_allsheets(S3.stART.sir.path)
S3.stART.hiv <- read_excel_allsheets(S3.stART.hiv.path)

S4.stART.inc <- read_excel_allsheets(S4.stART.inc.path)
S4.stART.ir <- read_excel_allsheets(S4.stART.ir.path)
S4.stART.sir <- read_excel_allsheets(S4.stART.sir.path)
S4.stART.hiv <- read_excel_allsheets(S4.stART.hiv.path)

colnames <- c("year", "median", "min", "max", "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9",
  "sim10", "sim11", "sim12", "sim13", "sim14", "sim15", "sim16", "sim17", "sim18", "sim19", "sim20",
  "sim21", "sim22", "sim23", "sim24", "sim25")

##########
# DEFINE DATA PROCESSING FUNCTIONS AND BASIC PLOT SETTINGS # ----
##########

combinesheets <- function(scenario){
  # Set column names
  scenario <- lapply(scenario, setNames, colnames)
  # Add outcome variables
  scenario[[1]]$outcome <- "CCC"
  scenario[[2]]$outcome <- "CCC"
  scenario[[3]]$outcome <- "CCC"
  scenario[[4]]$outcome <- "CCC"
  scenario[[5]]$outcome <- "CCC"
  scenario[[6]]$outcome <- "ICC"
  scenario[[7]]$outcome <- "ICC"
  scenario[[8]]$outcome <- "ICC"
  scenario[[9]]$outcome <- "ICC"
  scenario[[10]]$outcome <- "ICC"
  scenario[[11]]$outcome <- "OS"
  scenario[[12]]$outcome <- "OS"
  scenario[[13]]$outcome <- "OS"
  scenario[[14]]$outcome <- "OS"
  scenario[[15]]$outcome <- "OS"
  scenario[[16]]$outcome <- "SCRN"
  scenario[[17]]$outcome <- "SCRN"
  scenario[[18]]$outcome <- "SCRN"
  scenario[[19]]$outcome <- "SCRN"
  scenario[[20]]$outcome <- "SCRN"
 
  
  # Combine into one sheet
  out <- rbind.data.frame(scenario[[1]], scenario[[2]], scenario[[3]], scenario[[4]], scenario[[5]],
                          scenario[[6]], scenario[[7]], scenario[[8]], scenario[[9]], scenario[[10]],
                          scenario[[11]], scenario[[12]], scenario[[13]], scenario[[14]], scenario[[15]],
                          scenario[[16]], scenario[[17]], scenario[[18]], scenario[[19]], scenario[[20]])
  out[,c(2:6)] <- sapply(out[,c(2:6)], as.numeric)
  return(out)
}

# Function to round years, set column names, restrict to years 2000 to 2021, and combine estimates from all groups onto one sheet
dataprep.absnum <- function(scenario){
  # Set column names
  scenario <- lapply(scenario, setNames, colnames)
  # Restrict to years 2020 to 2121
  scenario <- lapply(scenario, function(x) filter(x, year >= 2000 & year <=2121))
  # Add group variables
  scenario[[1]]$group <- "Total"
  scenario[[2]]$group <- "HIVneg"
  scenario[[3]]$group <- "HIVpos_noART"
  scenario[[4]]$group <- "HIVpos_VS"
  scenario[[5]]$group <- "HIVpos_all"
  # Combine into one sheet
  out <- rbind.data.frame(scenario[[1]], scenario[[2]], scenario[[3]], scenario[[4]], scenario[[5]])
  out[,c(1:29)] <- sapply(out[,c(1:29)], as.numeric)
  
  return(out)
}

dataprep.hivprev <- function(scenario){
    # Set column names
    scenario <- lapply(scenario, setNames, colnames)
    # Restrict to years 2020 to 2121
    scenario <- lapply(scenario, function(x) filter(x, year >= 2000 & year <=2121))
    scenario[[1]]$group <- "Total" 
    out <- scenario[[1]]
    return(out)
}

# Function to calculate the percent reduction by group from 2020 over time
# pctreduc <- function(scenario){
#  pctreduc <- cbind.data.frame(year = scenario[,1], Total = scenario[,2]/scenario[1,2], HIVneg = scenario[,3]/scenario[1,3],
#                               HIVpos_noART = scenario[,4]/scenario[1,4], HIVpos_ART = scenario[,5]/scenario[1,5], HIVpos =  scenario[,6]/scenario[1,6]) 
#  return(pctreduc)
# }
# 
# Function to calculate the percent of cases attributable to WLHIV - for now set to use the median
calc_prophiv <- function(scenario){
  prophiv <- as.data.frame(matrix(ncol = 30, nrow = 122))
  colnames(prophiv) <- colnames(scenario)[1:29]
  prophiv$year <- seq(2000, 2121, 1)
  for(i in c(2000:2121)){
    prophiv[prophiv$year == i , c(2:29)] <- scenario[scenario$group == "HIVpos_all" & scenario$year == i, c(2:29)] / scenario[scenario$group == "Total" & scenario$year == i,  c(2:29)]
  }
  return(prophiv)
}

# Function to calculate the cumulative percent of cases attributable to WLHIV. prophiv.cuml = proportion of cases in WLHIV, propart.cuml = proportion of cases in virally suppressed WLHIV
calc_prophiv_cuml <- function(scenario, yearstart, yearend){
  scenario.trunc <- filter(scenario, year >=yearstart, year <= yearend)
  cuml <- scenario.trunc %>% group_by(group) %>% summarise(across(c(median:sim25), sum))
  prophiv.cuml <- cuml[cuml$group == "HIVpos_all", c(2:29)] / cuml[cuml$group == "Total", c(2:29)]
  prophiv.cuml$outcome = "Cuml proportion of cases in WLHIV"
  propart.cuml <- cuml[cuml$group == "HIVpos_VS", c(2:29)] / cuml[cuml$group == "Total", c(2:29)]
  propart.cuml$outcome = "Cuml proportion of cases in VS WLHIV"
  out <- rbind.data.frame(prophiv.cuml, propart.cuml)
  return(out)
}

## Plot settings
plot_background <- theme(panel.background = element_rect(fill="white", colour = "black")) + theme(panel.grid.major = element_line(colour = "grey90"))
plot_titles <- theme(plot.title = element_text(hjust = 0.5, size=14, colour = "black", face = "bold")) + 
  theme(legend.title = element_text(colour = "black", size = 13), legend.text = element_text(colour = "black", size = 13)) +
  theme(axis.title = element_text(colour = "black", size = 13), axis.text = element_text(colour = "black", size = 12))
colors <- RColorBrewer::brewer.pal(5, "Set2")


##########
# COMBINE ESTIMATES FOR EACH SCENARIO ONTO ONE SHEET # ----
##########

S0.stART.inc.comb <- dataprep.absnum(S0.stART.inc)
S0.stART.ir.comb <- dataprep.absnum(S0.stART.ir)
S0.stART.sir.comb <- dataprep.absnum(S0.stART.sir)
S0.stART.hiv.comb <- dataprep.hivprev(S0.stART.hiv)

S1.stART.inc.comb <- dataprep.absnum(S1.stART.inc)
S1.stART.ir.comb <- dataprep.absnum(S1.stART.ir)
S1.stART.sir.comb <- dataprep.absnum(S1.stART.sir)
S1.stART.hiv.comb <- dataprep.hivprev(S1.stART.hiv)

S2.stART.inc.comb <- dataprep.absnum(S2.stART.inc)
S2.stART.ir.comb <- dataprep.absnum(S2.stART.ir)
S2.stART.sir.comb <- dataprep.absnum(S2.stART.sir)
S2.stART.hiv.comb <- dataprep.hivprev(S2.stART.hiv)

S3.stART.inc.comb <- dataprep.absnum(S3.stART.inc)
S3.stART.ir.comb <- dataprep.absnum(S3.stART.ir)
S3.stART.sir.comb <- dataprep.absnum(S3.stART.sir)
S3.stART.hiv.comb <- dataprep.hivprev(S3.stART.hiv)

S4.stART.inc.comb <- dataprep.absnum(S4.stART.inc)
S4.stART.ir.comb <- dataprep.absnum(S4.stART.ir)
S4.stART.sir.comb <- dataprep.absnum(S4.stART.sir)
S4.stART.hiv.comb <- dataprep.hivprev(S4.stART.hiv)


##########
# CALCULATE IR in each year for each scenario # ----
##########
# S0.stART.ir.comb
# S0b.stART.ir.comb
# S1.stART.ir.comb
# S2.stART.ir.comb
# S3.stART.ir.comb

# S0.inART.ir.comb
# S0b.inART.ir.comb
# S1.inART.ir.comb
# S2.inART.ir.comb
# S3.inART.ir.comb
# 
# S0.stART.sir.comb
# S0b.stART.sir.comb
# S1.stART.sir.comb
# S2.stART.sir.comb
# S3.stART.sir.comb

# S0.inART.sir.comb
# S0b.inART.sir.comb
# S1.inART.sir.comb
# S2.inART.sir.comb
# S3.inART.sir.comb

# Combine into a DF for raw IR
S0.stART.ir.comb$scenario <- "Baseline cytology and vaccination"
S1.stART.ir.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.stART.ir.comb$scenario <- "S1 + more frequent screening for WLHIV"
S3.stART.ir.comb$scenario <- "S2 + 50% catch-up vaccination for WLHIV"
S4.stART.ir.comb$scenario <- "scenario description"

ir.comb.stART <- rbind.data.frame(S0.stART.ir.comb, S1.stART.ir.comb, S2.stART.ir.comb, S3.stART.ir.comb, S4.stART.ir.comb)
ir.comb.stART$scenario <- factor(ir.comb.stART$scenario, levels = c("Baseline cytology and vaccination", 
                                                        "Scaled up HPV testing and 90% vaccination",
                                                        "S1 + more frequent screening for WLHIV",
                                                        "S2 + 50% catch-up vaccination for WLHIV",
                                                        "scenario description"))
ir.comb.stART$outcome <- "Crude Incidence Rates"

# Combine into a DF for standardized IR
S0.stART.sir.comb$scenario <- "Baseline cytology and vaccination"
S1.stART.sir.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.stART.sir.comb$scenario <- "S1 + more frequent screening for WLHIV"
S3.stART.sir.comb$scenario <- "S2 + 50% catch-up vaccination for WLHIV"
S4.stART.sir.comb$scenario <- "scenario description"


sir.comb.stART <- rbind.data.frame(S0.stART.sir.comb, S1.stART.sir.comb, S2.stART.sir.comb, S3.stART.sir.comb, S4.stART.sir.comb)
sir.comb.stART$scenario <- factor(sir.comb.stART$scenario, levels = c("Baseline cytology and vaccination", 
                                                                    "Scaled up HPV testing and 90% vaccination",
                                                                    "S1 + more frequent screening for WLHIV",
                                                                    "S2 + 50% catch-up vaccination for WLHIV",
                                                                    "scenario description"))
sir.comb.stART$outcome <- "AS Incidence Rates"


##########
# Combine incident case counts for plotting # ----
##########

# Combine into a DF
S0.stART.inc.comb$scenario <- "Baseline cytology and vaccination"
S1.stART.inc.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.stART.inc.comb$scenario <- "S1 + more frequent screening for WLHIV"
S3.stART.inc.comb$scenario <- "S2 + 50% catch-up vaccination for WLHIV"
S4.stART.inc.comb$scenario <- "scenario description"

inc.comb.stART <- rbind.data.frame(S0.stART.inc.comb, S1.stART.inc.comb, S2.stART.inc.comb, S3.stART.inc.comb, S4.stART.inc.comb)
inc.comb.stART$scenario <- factor(inc.comb.stART$scenario, levels = c("Baseline cytology and vaccination", 
                                                        "Scaled up HPV testing and 90% vaccination",
                                                        "S1 + more frequent screening for WLHIV",
                                                        "S2 + 50% catch-up vaccination for WLHIV",
                                                        "scenario description"))
inc.comb.stART$outcome <- "Incident cases"


ratesandcases.stART <- rbind.data.frame(inc.comb.stART, ir.comb.stART, sir.comb.stART)

##########
# Combine HIV prevalence across scenarios # ----
##########

# Combine into a DF
S0.stART.hiv.comb$scenario <- "Baseline cytology and vaccination"
S1.stART.hiv.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.stART.hiv.comb$scenario <- "S1 + more frequent screening for WLHIV"
S3.stART.hiv.comb$scenario <- "S2 + 50% catch-up vaccination for WLHIV"
S4.stART.hiv.comb$scenario <- "scenario description"

hiv.comb.stART <- rbind.data.frame(S0.stART.hiv.comb, S1.stART.hiv.comb, S2.stART.hiv.comb, S3.stART.hiv.comb, S4.stART.hiv.comb)
hiv.comb.stART$scenario <- factor(hiv.comb.stART$scenario, levels = c("Baseline cytology and vaccination", 
                                                          "Scaled up HPV testing and 90% vaccination",
                                                          "S1 + more frequent screening for WLHIV",
                                                          "S2 + 50% catch-up vaccination for WLHIV",
                                                          "scenario description"))
hiv.comb.stART$outcome <- "HIV prevalence"



##########
# CALCULATE PAF in each year and cumulative PAF for each scenario # ----
##########

(prophiv.S0.stART <- calc_prophiv(S0.stART.inc.comb))
calc_prophiv_cuml(S0.stART.inc.comb, 2021, 2121)

(prophiv.S1.stART <- calc_prophiv(S1.stART.inc.comb))
calc_prophiv_cuml(S1.stART.inc.comb, 2021, 2121)

(prophiv.S2.stART <- calc_prophiv(S2.stART.inc.comb))
calc_prophiv_cuml(S2.stART.inc.comb, 2021, 2121)

(prophiv.S3.stART <- calc_prophiv(S3.stART.inc.comb))
calc_prophiv_cuml(S3.stART.inc.comb, 2021, 2121)

(prophiv.S4.stART <- calc_prophiv(S4.stART.inc.comb))
calc_prophiv_cuml(S4.stART.inc.comb, 2021, 2121)

# Combine into a data frame
prophiv.S0.stART$scenario <- "Baseline cytology and vaccination"
prophiv.S1.stART$scenario <- "Scaled up HPV testing and 90% vaccination"
prophiv.S2.stART$scenario <- "S1 + more frequent screening for WLHIV"
prophiv.S3.stART$scenario <- "S2 + 50% catch-up vaccination for WLHIV"
prophiv.S4.stART$scenario <- "scenario description"

prophiv.stART.comb <- rbind.data.frame(prophiv.S0.stART, prophiv.S1.stART, prophiv.S2.stART, prophiv.S3.stART, prophiv.S4.stART)
prophiv.stART.comb$scenario <- factor(prophiv.stART.comb$scenario, levels = c("Baseline cytology and vaccination", 
                                                          "Scaled up HPV testing and 90% vaccination",
                                                          "S1 + more frequent screening for WLHIV",
                                                          "S2 + 50% catch-up vaccination for WLHIV",
                                                          "scenario description"))


############################
## PLOTS
############################

colors <- RColorBrewer::brewer.pal(5, "Set2")

## Look at *median* outcomes across scenarios.
sir.comb.stART %>%  # Change the name of the dataframe to different outcomes (ir.comb.stART, sir.comb.stART, inc.comb.stART, )
  filter(group %in% c("Total", "HIVneg", "HIVpos_noART", "HIVpos_VS")) %>%
  ggplot() + 
    geom_line(aes(x=year, y = median, colour = scenario)) +
    facet_wrap(vars(group), scales = "free")

ggplot(prophiv.stART.comb) + 
    geom_line(aes(x=year, y = median, colour = scenario)) +
    #geom_line(data = hiv.comb.inART, aes(x=year, y = median, colour = scenario)) +
    geom_line(data = hiv.comb.stART, aes(x=year, y = median, colour = scenario), linetype = "dashed")


### Incidence rates
## Building up the scenarios for presentation

# Scenario 0 - Total pop only
#png("irS0_tot.png", width = 800, height = 400)
pdf("irS0_tot.pdf", width = 10, height = 5)
S0.stART.ir.comb %>%
  filter(group == "Total") %>%
  ggplot() + 
    geom_line(aes(x = year, y = median), colour = "black", size = 1.25) +
    scale_y_continuous(breaks = seq(0, 105, 25), limits = c(0, 115)) +
    labs(x = "year", y = "Incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black")) +
    plot_background 
dev.off()

# Scenario 0 - Adding in HIV states
#png("irS0_byhiv.png", width = 800, height = 400)
pdf("irS0_byhiv.pdf", width = 10, height = 5) # legend position for this size shoudl be c(0.85, 0.85)
#pdf("irS0_byhiv.small.pdf", width = 6, height = 5)  # save it again at half size, and change legend position to c(0.8, 0.85)
S0.stART.ir.comb %>%
  filter(group %in% c("Total", "HIVneg", "HIVpos_all")) %>%
  ggplot() + 
    geom_line(aes(x = year, y = median, colour = group), size = 1.25) +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos_all"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
    labs(x = "year", y = "Incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()


===> STOPPED UPDATING HERE


# png("irS0_byhivandart.png", width = 800, height = 400)
# ggplot(S0b.stART.ir.comb) + 
#     geom_line(aes(x = year, y = Total, colour = "Total"), size = 1.25) +
#     geom_line(aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1) +
#     geom_line(aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1) +
#     geom_line(aes(x = year, y = HIVpos_noART, colour = "HIVpos_noART"), size = 1) +
#     geom_line(aes(x = year, y = HIVpos_ART, colour = "HIVpos_ART"), size = 1) +
#     scale_colour_manual(values = c("Total" = "black", "HIVneg" = colors[5], "HIVpos" = colors[3], "HIVpos_noART" = colors[2], "HIVpos_ART" = colors[4]),
#                         breaks = c("Total", "HIVneg", "HIVpos", "HIVpos_noART", "HIVpos_ART"),
#                         labels = c("All females", "HIV-negative", "HIV-positive (all)", "HIV-positive, untreated", "HIV-positive, virally suppressed"),
#                         name = NULL) +
#     labs(x = "year", y = "Incidence rate per 100,000 women") +
#     theme(axis.title = element_text(size = 18),
#           axis.text = element_text(size = 16),
#           axis.line = element_line(color = "black"),
#           legend.text = element_text(size = 18),
#           legend.position = c(0.85, 0.85)) +
#     plot_background
# dev.off()

# Add in scenario with increasing ART
#png("irS0b_byhiv_inART.png", width = 800, height = 400)
pdf("irS0b_byhiv_inART.pdf", width = 6, height = 5)
ggplot() + 
    geom_line(data = S0b.stART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.75, alpha = 0.55) +
    geom_line(data = S0b.stART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1.5, alpha = 0.55) +
    geom_line(data = S0b.stART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1.5, alpha = 0.55) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1) +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 105, 25), limits = c(0, 115)) +
    labs(x = "year", y = "Incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()

# Add in S1 with increasing ART
#png("irS0bS1_byhiv.png", width = 800, height = 400)
pdf("irS0bS1_byhiv.pdf", width = 6, height = 5)
ggplot() + 
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1) +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25, linetype = "dashed") +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1, linetype = "dashed") +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1, linetype = "dashed") +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 105, 25), limits = c(0, 115)) +
    labs(x = "year", y = "Incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()


# Add in S3 with increasing ART
#png("irS0bS1S3_byhiv.png", width = 800, height = 400)
#pdf("irS0bS1S3_byhiv.pdf", width = 10, height = 5)
pdf("irS0bS1S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1) +
    geom_line(data = S0b.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1) +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25, linetype = "dashed") +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1, linetype = "dashed") +
    geom_line(data = S1.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1, linetype = "dashed") +
    geom_line(data = S3.inART.ir.comb, aes(x = year, y = Total, colour = "Total"), size = 1.25, linetype = "dotdash") +
    geom_line(data = S3.inART.ir.comb, aes(x = year, y = HIVneg, colour = "HIVneg"), size = 1, linetype = "dotdash") +
    geom_line(data = S3.inART.ir.comb, aes(x = year, y = HIVpos, colour = "HIVpos"), size = 1, linetype = "dotdash") +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 105, 25), limits = c(0, 115)) +
    labs(x = "year", y = "Incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()


### PAF and HIV prevalence

## Scenario S0b
pdf("paf.s0b.small.pdf", width = 6, height = 5) 
ggplot(paf.S0b.stART) +
    geom_area(aes(x = year, y = paf_hiv), fill = colors[3]) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(x = "year", y = "Proportion of cases among WLHIV") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

# Add HIV prevalence
pdf("paf.prev.s0b.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = paf.S0b.stART, aes(x = year, y = paf_hiv), fill = colors[3]) +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") + 
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

## Scenario S0b with increasing ART
pdf("paf.prev.s0b.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = paf.S0b.stART, aes(x = year, y = paf_hiv), fill = colors[3], alpha = 0.6) +
    geom_area(data = paf.S0b.inART, aes(x = year, y = paf_hiv), fill = colors[3], alpha = 0.6) +
    geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#486199") +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = Prevalence), fill = "#8B7DCA") + 
    geom_area(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") + 
    geom_line(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), colour = "#3C2F74") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

#Highlight the initial increase in PAF and then eventual decrease
S0b.min <- pmin(paf.S0b.stART$paf_hiv, paf.S0b.inART$paf_hiv)
paf.S0b.min <- cbind.data.frame(year = paf.S0b.stART$year, min_paf = S0b.min)

pdf("paf.prev.s0b.inART.showincr.small.pdf", width = 6, height = 5) 
ggplot() +
    geom_area(data = paf.S0b.inART, aes(x = year, y = paf_hiv), fill = "#A4F283") +
    geom_area(data = paf.S0b.stART, aes(x = year, y = paf_hiv), fill = colors[3], alpha = 0.6) +
    geom_area(data = paf.S0b.min, aes(x = year, y = min_paf), fill = colors[3]) +
    geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#486199") +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = Prevalence), fill = "#8B7DCA") + 
    geom_area(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") + 
    geom_line(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), colour = "#3C2F74") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

pdf("paf.prev.s0b.inART.showdecr.small.pdf", width = 6, height = 5) 
ggplot() +
    geom_area(data = paf.S0b.inART, aes(x = year, y = paf_hiv), fill = "#A4F283") +
    geom_area(data = paf.S0b.stART, aes(x = year, y = paf_hiv), fill = "#C76F98", alpha = 0.8) +
    geom_area(data = paf.S0b.min, aes(x = year, y = min_paf), fill = colors[3]) +
    geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#486199") +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = Prevalence), fill = "#8B7DCA") + 
    geom_area(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") + 
    geom_line(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), colour = "#3C2F74") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

## Add in S1
pdf("paf.prev.s0bs1.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = paf.S0b.inART, aes(x = year, y = paf_hiv, fill = "Baseline")) +
    geom_area(data = paf.S1.inART, aes(x = year, y = paf_hiv, fill = "S1"), alpha = 0.55) +
    geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#486199") +
    geom_line(data = paf.S1.inART, aes(x = year, y = paf_hiv), colour = "#486199", linetype = "dashed") +
    # geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#738ABF") +
    # geom_line(data = paf.S1.inART, aes(x = year, y = paf_hiv), colour = "#7AB4F5") +
    geom_area(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_fill_manual(values = c("Baseline" = colors[3], "S1" = "#97C4F7"),
                                          breaks = c("Baseline", "S1"),
                                          labels = c("Baseline screening and vaccination", "Scaled up HPV testing and 90% vax"),
                                          name = NULL) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = "none") +
    plot_background
dev.off()

## Add in S3
pdf("paf.prev.s0bs1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = paf.S0b.inART, aes(x = year, y = paf_hiv, fill = "Baseline")) +
    geom_area(data = paf.S1.inART, aes(x = year, y = paf_hiv, fill = "S1"), alpha = 0.55) +
    geom_area(data = paf.S3.inART, aes(x = year, y = paf_hiv, fill = "S3")) +
    geom_line(data = paf.S0b.inART, aes(x = year, y = paf_hiv), colour = "#486199") +
    geom_line(data = paf.S1.inART, aes(x = year, y = paf_hiv), colour = "#486199", linetype = "dashed") +
    geom_line(data = paf.S3.inART, aes(x = year, y = paf_hiv), colour = "#486199", linetype = "dotdash") +
    geom_area(data = S0b.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#5F4BB6") +
    # geom_area(data = S1.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#A15AD1") + 
    # geom_area(data = S3.inART.hiv.comb, aes(x = year, y = Prevalence), fill = "#8771EB") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_fill_manual(values = c("Baseline" = colors[3], "S1" = "#97C4F7", "S3" = "#A0C9F8"),
                      breaks = c("Baseline", "S1", "S3"),
                      labels = c("Baseline screening and vaccination", "Scaled up HPV testing and 90% vax",
                                 "Increased HPV screening and catch-up vax for WLHIV"),
                      name = NULL) +
    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 18),
          legend.position = "none") +
    plot_background
dev.off()

# ggplot(paf.inART.comb) + 
#     geom_line(aes(x = year, y = paf_hiv, colour = scenario))
# ggplot(hiv.comb.inART) +
#     geom_line(aes(x = year, y = Prevalence, colour = scenario))

    
## Extract statistics
# Incidence rates in 2000, 2020, and 2070 and % change
(startend.stART <- filter(S0b.start.ir.comb, year %in% c(2000, 2020, 2070), scenario == "Baseline w/ 2x cytology"))
startend.stART %>% group_by(scenario) %>% mutate(pctchange_gen_2020 = Total[year == 2020] / Total[year == 2000],
                                                 pctchange_gen_2070 = Total[year == 2070] / Total[year == 2020],
                                                 pctchange_neg_2020 = HIVneg[year == 2020] / HIVneg[year == 2000],
                                                 pctchange_neg_2070 = HIVneg[year == 2070] / HIVneg[year == 2020],
                                                 pctchange_pos_2020 = HIVpos[year == 2020] / HIVpos[year == 2000],
                                                 pctchange_pos_2070 = HIVpos[year == 2070] / HIVpos[year == 2020]) %>%
    filter(year == 2000) %>%
    select(scenario, pctchange_gen_2020, pctchange_gen_2070, pctchange_neg_2020, pctchange_neg_2070, pctchange_pos_2020, pctchange_pos_2070)
    

(startend.inART <- filter(ir.comb.inART, year %in% c(2000, 2020, 2070), scenario %in% c("Baseline w/ 2x cytology", 
                                                                                 "Scaled up HPV testing and 90% vaccination",
                                                                                 "S2 + 50% catch-up vaccination for WLHIV")))
startend.inART %>% group_by(scenario) %>% mutate(pctchange_gen_2020 = Total[year == 2020] / Total[year == 2000],
                                                 pctchange_gen_2070 = Total[year == 2070] / Total[year == 2020],
                                           pctchange_neg_2020 = HIVneg[year == 2020] / HIVneg[year == 2000],
                                           pctchange_neg_2070 = HIVneg[year == 2070] / HIVneg[year == 2020],
                                           pctchange_pos_2020 = HIVpos[year == 2020] / HIVpos[year == 2000],
                                           pctchange_pos_2070 = HIVpos[year == 2070] / HIVpos[year == 2020]) %>%
    filter(year == 2000) %>%
    select(scenario, pctchange_gen_2020, pctchange_gen_2070, pctchange_neg_2020, pctchange_neg_2070, pctchange_pos_2020, pctchange_pos_2070)
   
    # Pctchange relative to baseline/other scenarios
        (pctreducS1S0b <- cbind.data.frame(Total = 1 - (startend.inART$Total[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070] / 
                                           startend.inART$Total[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070]),
                                        HIVneg = 1 - (startend.inART$HIVneg[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070] / 
                                           startend.inART$HIVneg[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070]),
                                        HIVpos = 1 - (startend.inART$HIVpos[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070] / 
                                            startend.inART$HIVpos[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070])))
        (pctreducS3S0b <- cbind.data.frame(Total = 1 - (startend.inART$Total[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$Total[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070]),
                                     HIVneg = 1 - (startend.inART$HIVneg[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$HIVneg[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070]),
                                    HIVpos = 1 - (startend.inART$HIVpos[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$HIVpos[startend.inART$scenario == "Baseline w/ 2x cytology" & startend.inART$year == 2070])))
        (pctreducS3S1 <- cbind.data.frame(Total = 1 - (startend.inART$Total[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$Total[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070]),
                                    HIVneg = 1 - (startend.inART$HIVneg[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$HIVneg[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070]),
                                    HIVpos = 1 - (startend.inART$HIVpos[startend.inART$scenario == "S2 + 50% catch-up vaccination for WLHIV" & startend.inART$year == 2070] / 
                                        startend.inART$HIVpos[startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & startend.inART$year == 2070])))
        
    # find the peak year for each group
    (max.inc <- ir.comb.stART %>% filter(scenario == "Baseline w/ 2x cytology") %>% summarise_at(vars(Total, HIVneg, HIVpos), max))
    ir.comb.stART$year[ir.comb.stART$Total == max.inc$Total]
    ir.comb.stART$year[ir.comb.stART$HIVneg == max.inc$HIVneg]
    ir.comb.stART$year[ir.comb.stART$HIVpos == max.inc$HIVpos]
    
    (max.inc <- ir.comb.inART %>% filter(scenario == "Baseline w/ 2x cytology") %>% summarise_at(vars(Total, HIVneg, HIVpos), max))
    ir.comb.inART$year[ir.comb.inART$Total == max.inc$Total]
    ir.comb.inART$year[ir.comb.inART$HIVneg == max.inc$HIVneg]
    ir.comb.inART$year[ir.comb.inART$HIVpos == max.inc$HIVpos]

    # To find how much higher incidence is in the ART scale-up scenario than the sceanrio without scale-up. Adjust group, scenario, and year as needed
    ir.comb.inART$HIVpos[ir.comb.inART$scenario == "Baseline w/ 2x cytology" & ir.comb.inART$year==2070] / 
        ir.comb.stART$HIVpos[ir.comb.stART$scenario== "Baseline w/ 2x cytology" & ir.comb.stART$year==2070]

# PAF in 2000, 2020, and 2070 and % change
(paf.startend.stART <- filter(paf.stART.comb, year %in% c(2000, 2020, 2070), scenario == "Baseline w/ 2x cytology") %>% select(year, scenario, paf_hiv))
paf.startend.stART %>% group_by(scenario) %>% mutate(pctchange_2020 = paf_hiv[year == 2020] / paf_hiv[year == 2000],
                                                 pctchange_2070 = paf_hiv[year == 2070] / paf_hiv[year == 2020]) %>%
    filter(year == 2000) %>%
    select(scenario, pctchange_2020, pctchange_2070)

(paf.startend.inART <- filter(paf.inART.comb, year %in% c(2000, 2020, 2070), scenario %in% c("Baseline w/ 2x cytology", 
                                                                                        "Scaled up HPV testing and 90% vaccination",
                                                                                        "S2 + 50% catch-up vaccination for WLHIV")) %>% 
                    select(year, scenario, paf_hiv))
paf.startend.inART %>% group_by(scenario) %>% mutate(pctchange_2020 = paf_hiv[year == 2020] / paf_hiv[year == 2000],
                                                     pctchange_2070 = paf_hiv[year == 2070] / paf_hiv[year == 2020]) %>%
    filter(year == 2000) %>%
    select(scenario, pctchange_2020, pctchange_2070)

    # find the peak year for each group
    (max.inc <- paf.stART.comb %>% filter(scenario == "Baseline w/ 2x cytology") %>% summarise_at(vars(paf_hiv), max))
    paf.stART.comb$year[paf.stART.comb$paf_hiv == max.inc$paf_hiv]
 
    (max.inc <- paf.inART.comb %>% group_by(scenario) %>% summarise_at(vars(paf_hiv), max))[-3 ,]
    paf.inART.comb$year[paf.inART.comb$scenario == "Baseline w/ 2x cytology" & paf.inART.comb$paf_hiv == max.inc$paf_hiv[1]]
    paf.inART.comb$year[paf.inART.comb$scenario == "Scaled up HPV testing and 90% vaccination" & paf.inART.comb$paf_hiv == max.inc$paf_hiv[2]]
    paf.inART.comb$year[paf.inART.comb$scenario == "S2 + 50% catch-up vaccination for WLHIV" & paf.inART.comb$paf_hiv == max.inc$paf_hiv[4]]
    
# HIV prevalence in 2000, 2020, and 2070 and % change
(hiv.startend.stART <- filter(hiv.comb.stART, year %in% c(2000, 2020, 2070), scenario == "Baseline w/ 2x cytology") %>% select(year, scenario, Prevalence))
hiv.startend.stART %>% group_by(scenario) %>% mutate(pctchange_2020 = Prevalence[year == 2020] / Prevalence[year == 2000],
                                                    pctchange_2070 = Prevalence[year == 2070] / Prevalence[year == 2020]) %>%
        filter(year == 2000) %>%
        select(scenario, pctchange_2020, pctchange_2070)
    
(hiv.startend.inART <- filter(hiv.comb.inART, year %in% c(2000, 2020, 2070), scenario %in% c("Baseline w/ 2x cytology", 
                                                                                            "Scaled up HPV testing and 90% vaccination",
                                                                                            "S2 + 50% catch-up vaccination for WLHIV")) %>% 
            select(year, scenario, Prevalence))
hiv.startend.inART %>% group_by(scenario) %>% mutate(pctchange_2020 = Prevalence[year == 2020] / Prevalence[year == 2000],
                                                     pctchange_2070 = Prevalence[year == 2070] / Prevalence[year == 2020]) %>%
        filter(year == 2000) %>%
        select(scenario, pctchange_2020, pctchange_2070)
    
# find the peak year for each group
(max.inc <- hiv.comb.stART %>% filter(scenario == "Baseline w/ 2x cytology") %>% summarise_at(vars(Prevalence), max))
hiv.comb.stART$year[hiv.comb.stART$Prevalence == max.inc$Prevalence]
    
(max.inc <- hiv.comb.inART %>% filter(scenario == "Baseline w/ 2x cytology") %>% summarise_at(vars(Prevalence), max))
hiv.comb.inART$year[hiv.comb.inART$Prevalence == max.inc$Prevalence]


# Disparities over time
startend.stART %>% group_by(scenario) %>% mutate(disparity_2000 = HIVpos[year == 2000] / HIVneg[year == 2000],
                                                 disparity_2020 = HIVpos[year == 2020] / HIVneg[year == 2020],
                                                 disparity_2070 = HIVpos[year == 2070] / HIVneg[year == 2070]) %>%
    filter(year == 2000) %>%
    select(scenario, disparity_2000, disparity_2020, disparity_2070)

startend.inART %>% group_by(scenario) %>% mutate(disparity_2000 = HIVpos[year == 2000] / HIVneg[year == 2000],
                                                 disparity_2020 = HIVpos[year == 2020] / HIVneg[year == 2020],
                                                 disparity_2070 = HIVpos[year == 2070] / HIVneg[year == 2070]) %>%
    filter(year == 2000) %>%
    select(scenario, disparity_2000, disparity_2020, disparity_2070)

# Disparities in terms of PAF_HIV divided by HIV prev
disparity.stART <- cbind.data.frame(year = paf.stART.comb$year, scenario = paf.stART.comb$scenario, paf = paf.stART.comb$paf_hiv, prevalence = hiv.comb.stART$Prevalence)
disparity.inART <- cbind.data.frame(year = paf.inART.comb$year, scenario = paf.inART.comb$scenario, paf = paf.inART.comb$paf_hiv, prevalence = hiv.comb.inART$Prevalence)

disparity.stART %>% group_by(scenario) %>% mutate(disparity_2000 = paf[year == 2000] / prevalence[year == 2000],
                                                 disparity_2020 = paf[year == 2020] / prevalence[year == 2020],
                                                 disparity_2070 = paf[year == 2070] / prevalence[year == 2070]) %>%
    filter(year == 2000) %>%
    select(scenario, disparity_2000, disparity_2020, disparity_2070)

disparity.inART %>% group_by(scenario) %>% mutate(disparity_2000 = paf[year == 2000] / prevalence[year == 2000],
                                                  disparity_2020 = paf[year == 2020] / prevalence[year == 2020],
                                                  disparity_2070 = paf[year == 2070] / prevalence[year == 2070]) %>%
    filter(year == 2000) %>%
    select(scenario, disparity_2000, disparity_2020, disparity_2070)
