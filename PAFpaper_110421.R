### Script to process output on one-time screening optimal age ###

# This script takes output from the simulation model and processes it to produce plots and statistics for tables and text. 
# The script expects four output files from each scenario:
# 1. Cervical cancer cases
# 2. Cervical cancer incidence rates (raw)
# 3. Cervical cancer incidence rates (age-standardized)
# 4. HIV prevalence in females 
# For each, there should be separate tabs/sheets for each of the following groups: the Total population, HIV-negative women, HIV-positive untreated women, HIV-positive women on ART, all HIv-positive women __IN THAT ORDER__. 
# Within each tab column 1 should show the year of simulation from 2020 to 2070. Column 2 should show the outcome for the relevant year, group, and scenario

## The filepaths and names correspond to these scenarios:
# S0b = Cytology at age 35 with 48% coverage throughout the timeframe. 9v HPV vaccination at 57%. No ART scale-up.
# S0 = Cytology at age 35 with 48% coverage throughout the timeframe. 9v HPV vaccination at 57%. ART scale-up to 90-90-90 by 2030.
# S1 = HPV screening at ages 35 and 45 scaled up to 90% by 2045 throughout the timeframe. 9v HPV vaccination at 90%. ART scale-up to 90-90-90 by 2030.
# S2 = S1 but with catch-up vaccination for WLHIV ages 15-24 at 50%. ART scale-up to 90-90-90 by 2030.
# S3 = S2 but with screening for WLHIV every 3 years ages 25-49. ART scale-up to 90-90-90 by 2030.

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
setwd("H:/HHCoM/HHCoM_Results/22Apr20Ph2V11_2v57BaseVax_spCytoScreen_noVMMChpv_hivInt2017_pHIV-S0_6_1")

## Set filepaths <<--
# List file names for files containing incidence output for each scenario
# These files are expected to have one tab for each population subgroup (All, HIV-, untreated HIV+, treated HIV+, all HIV+)
# and for each tab to list the years in column 1 and the outcome in column 2
S0b.stART.inc.path <- "PAF_crudeAnnualCC_S0.xlsx"
S0b.stART.ir.path <- "PAF_crudeICC_S0.xlsx"
S0b.stART.sir.path <- "PAF_ASICC_S0.xlsx"
S0b.stART.hiv.path <- "PAF_crudeHivFaged15plus_S0.xlsx"

S0.inART.inc.path <- "PAF_crudeAnnualCC_S1.xlsx"
S0.inART.ir.path <- "PAF_crudeICC_S1.xlsx"
S0.inART.sir.path <- "PAF_ASICC_S1.xlsx"
S0.inART.hiv.path <- "PAF_crudeHivFaged15plus_S1.xlsx"

S1.inART.inc.path <- "PAF_crudeAnnualCC_S2.xlsx"
S1.inART.ir.path <- "PAF_crudeICC_S2.xlsx"
S1.inART.sir.path <- "PAF_ASICC_S2.xlsx"
S1.inART.hiv.path <- "PAF_crudeHivFaged15plus_S2.xlsx"

S2.inART.inc.path <- "PAF_crudeAnnualCC_S3.xlsx"
S2.inART.ir.path <- "PAF_crudeICC_S3.xlsx"
S2.inART.sir.path <- "PAF_ASICC_S3.xlsx"
S2.inART.hiv.path <- "PAF_crudeHivFaged15plus_S3.xlsx"

S3.inART.inc.path <- "PAF_crudeAnnualCC_S4.xlsx"
S3.inART.ir.path <- "PAF_crudeICC_S4.xlsx"
S3.inART.sir.path <- "PAF_ASICC_S4.xlsx"
S3.inART.hiv.path <- "PAF_crudeHivFaged15plus_S4.xlsx"

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

S0b.stART.inc <- read_excel_allsheets(S0b.stART.inc.path)
S0b.stART.ir <- read_excel_allsheets(S0b.stART.ir.path)
S0b.stART.sir <- read_excel_allsheets(S0b.stART.sir.path)
S0b.stART.hiv <- read_excel_allsheets(S0b.stART.hiv.path)

S0.inART.inc <- read_excel_allsheets(S0.inART.inc.path)
S0.inART.ir <- read_excel_allsheets(S0.inART.ir.path)
S0.inART.sir <- read_excel_allsheets(S0.inART.sir.path)
S0.inART.hiv <- read_excel_allsheets(S0.inART.hiv.path)

S1.inART.inc <- read_excel_allsheets(S1.inART.inc.path)
S1.inART.ir <- read_excel_allsheets(S1.inART.ir.path)
S1.inART.sir <- read_excel_allsheets(S1.inART.sir.path)
S1.inART.hiv <- read_excel_allsheets(S1.inART.hiv.path)

S2.inART.inc <- read_excel_allsheets(S2.inART.inc.path)
S2.inART.ir <- read_excel_allsheets(S2.inART.ir.path)
S2.inART.sir <- read_excel_allsheets(S2.inART.sir.path)
S2.inART.hiv <- read_excel_allsheets(S2.inART.hiv.path)

S3.inART.inc <- read_excel_allsheets(S3.inART.inc.path)
S3.inART.ir <- read_excel_allsheets(S3.inART.ir.path)
S3.inART.sir <- read_excel_allsheets(S3.inART.sir.path)
S3.inART.hiv <- read_excel_allsheets(S3.inART.hiv.path)

colnames <- c("year", "median", "min", "max", "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9",
  "sim10", "sim11", "sim12", "sim13", "sim14", "sim15", "sim16", "sim17", "sim18", "sim19", "sim20",
  "sim21", "sim22", "sim23", "sim24", "sim25")

##########
# DEFINE DATA PROCESSING FUNCTIONS AND BASIC PLOT SETTINGS # ----
##########

#combinesheets <- function(scenario){
#  # Set column names
#  scenario <- lapply(scenario, setNames, colnames)
#  # Add outcome variables
#  scenario[[1]]$outcome <- "CCC"
#  scenario[[2]]$outcome <- "CCC"
#  scenario[[3]]$outcome <- "CCC"
#  scenario[[4]]$outcome <- "CCC"
#  scenario[[5]]$outcome <- "CCC"
#  scenario[[6]]$outcome <- "ICC"
#  scenario[[7]]$outcome <- "ICC"
#  scenario[[8]]$outcome <- "ICC"
#  scenario[[9]]$outcome <- "ICC"
#  scenario[[10]]$outcome <- "ICC"
#  scenario[[11]]$outcome <- "OS"
#  scenario[[12]]$outcome <- "OS"
#  scenario[[13]]$outcome <- "OS"
#  scenario[[14]]$outcome <- "OS"
#  scenario[[15]]$outcome <- "OS"
#  scenario[[16]]$outcome <- "SCRN"
#  scenario[[17]]$outcome <- "SCRN"
#  scenario[[18]]$outcome <- "SCRN"
#  scenario[[19]]$outcome <- "SCRN"
#  scenario[[20]]$outcome <- "SCRN"
# 
#  
#  # Combine into one sheet
#  out <- rbind.data.frame(scenario[[1]], scenario[[2]], scenario[[3]], scenario[[4]], scenario[[5]],
#                          scenario[[6]], scenario[[7]], scenario[[8]], scenario[[9]], scenario[[10]],
#                          scenario[[11]], scenario[[12]], scenario[[13]], scenario[[14]], scenario[[15]],
#                          scenario[[16]], scenario[[17]], scenario[[18]], scenario[[19]], scenario[[20]])
#  out[,c(2:6)] <- sapply(out[,c(2:6)], as.numeric)
#  return(out)
#}

# Function to round years, set column names, restrict to years 2000 to 2121, and combine estimates from all groups onto one sheet
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
 
# Function to calculate the proportion of cases among WLHIV
calc_prophiv <- function(scenario){
  prophiv <- as.data.frame(matrix(ncol = 29, nrow = 122))
  colnames(prophiv) <- colnames(scenario)[1:29]
  prophiv$year <- seq(2000, 2121, 1)
  for(i in c(2000:2121)){
    prophiv[prophiv$year == i , c(5:29)] <- scenario[scenario$group == "HIVpos_all" & scenario$year == i, c(5:29)] / scenario[scenario$group == "Total" & scenario$year == i,  c(5:29)]
    prophiv[prophiv$year == i , c(2:4)] <- c(apply(prophiv[prophiv$year == i , c(5:29)] , 1 , median) , 
                                             apply(prophiv[prophiv$year == i , c(5:29)] , 1 , min) , 
                                             apply(prophiv[prophiv$year == i , c(5:29)] , 1 , max))
  }
  return(prophiv)
}

# Function to calculate the proportion of cases among WLHIV divided by HIV prevalence
calc_ccHivRatio <- function(propCasesHiv , hivPrev){
  ccHivRatio <- as.data.frame(matrix(ncol = 29, nrow = 122))
  colnames(ccHivRatio) <- colnames(propCasesHiv)[1:29]
  ccHivRatio$year <- seq(2000, 2121, 1)
  for(i in c(2000:2121)){
    ccHivRatio[ccHivRatio$year == i , c(5:29)] <- propCasesHiv[propCasesHiv$year == i, c(5:29)] / hivPrev[hivPrev$year == i,  c(5:29)]
    ccHivRatio[ccHivRatio$year == i , c(2:4)] <- c(apply(ccHivRatio[ccHivRatio$year == i , c(5:29)] , 1 , median) , 
                                             apply(ccHivRatio[ccHivRatio$year == i , c(5:29)] , 1 , min) , 
                                             apply(ccHivRatio[ccHivRatio$year == i , c(5:29)] , 1 , max))
  }
  return(ccHivRatio)
}

# Function to calculate the percent of cases among WLHIV with VS
calc_prophivVS <- function(scenario){
  prophivVS <- as.data.frame(matrix(ncol = 29, nrow = 122))
  colnames(prophivVS) <- colnames(scenario)[1:29]
  prophivVS$year <- seq(2000, 2121, 1)
  for(i in c(2000:2121)){
    prophivVS[prophivVS$year == i , c(5:29)] <- scenario[scenario$group == "HIVpos_VS" & scenario$year == i, c(5:29)] / scenario[scenario$group == "Total" & scenario$year == i,  c(5:29)]
    prophivVS[prophivVS$year == i , c(2:4)] <- c(apply(prophivVS[prophivVS$year == i , c(5:29)] , 1 , median) , 
                                             apply(prophivVS[prophivVS$year == i , c(5:29)] , 1 , min) , 
                                             apply(prophivVS[prophivVS$year == i , c(5:29)] , 1 , max))
  }
  return(prophivVS)
}

# Function to calculate the percent of cases among WLHIV among WLHIV with VS
calc_prophivVSofHIV <- function(scenario){
  prophivVSofHIV <- as.data.frame(matrix(ncol = 29, nrow = 122))
  colnames(prophivVSofHIV) <- colnames(scenario)[1:29]
  prophivVSofHIV$year <- seq(2000, 2121, 1)
  for(i in c(2000:2121)){
    prophivVSofHIV[prophivVSofHIV$year == i , c(5:29)] <- scenario[scenario$group == "HIVpos_VS" & scenario$year == i, c(5:29)] / scenario[scenario$group == "HIVpos_all" & scenario$year == i,  c(5:29)]
    prophivVSofHIV[prophivVSofHIV$year == i , c(2:4)] <- c(apply(prophivVSofHIV[prophivVSofHIV$year == i , c(5:29)] , 1 , median) , 
                                                 apply(prophivVSofHIV[prophivVSofHIV$year == i , c(5:29)] , 1 , min) , 
                                                 apply(prophivVSofHIV[prophivVSofHIV$year == i , c(5:29)] , 1 , max))
  }
  return(prophivVSofHIV)
}

# Function to calculate the cumulative proportion of cases among WLHIV. prophiv.cuml = proportion of cases in WLHIV, propart.cuml = proportion of cases in virally suppressed WLHIV
calc_prophiv_cuml <- function(scenario, yearstart, yearend){
  scenario.trunc <- filter(scenario, year >=yearstart, year <= yearend)
  cuml <- scenario.trunc %>% group_by(group) %>% summarise(across(c(median:sim25), sum))
  prophiv.cuml <- cuml[cuml$group == "HIVpos_all", c(2:29)] / cuml[cuml$group == "Total", c(2:29)]
  prophiv.cuml[1,c(1:3)] <- c(apply(prophiv.cuml[1,c(4:28)] , 1 , median) , 
                               apply(prophiv.cuml[1,c(4:28)] , 1 , min) , 
                               apply(prophiv.cuml[1,c(4:28)] , 1 , max))
  prophiv.cuml$outcome = "Cuml proportion of cases in WLHIV"
  propart.cuml <- cuml[cuml$group == "HIVpos_VS", c(2:29)] / cuml[cuml$group == "Total", c(2:29)]
  propart.cuml[1,c(1:3)] <- c(apply(propart.cuml[1,c(4:28)] , 1 , median) , 
                            apply(propart.cuml[1,c(4:28)] , 1 , min) , 
                            apply(propart.cuml[1,c(4:28)] , 1 , max))
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

S0b.stART.inc.comb <- dataprep.absnum(S0b.stART.inc)
S0b.stART.ir.comb <- dataprep.absnum(S0b.stART.ir)
S0b.stART.sir.comb <- dataprep.absnum(S0b.stART.sir)
S0b.stART.hiv.comb <- dataprep.hivprev(S0b.stART.hiv)

S0.inART.inc.comb <- dataprep.absnum(S0.inART.inc)
S0.inART.ir.comb <- dataprep.absnum(S0.inART.ir)
S0.inART.sir.comb <- dataprep.absnum(S0.inART.sir)
S0.inART.hiv.comb <- dataprep.hivprev(S0.inART.hiv)

S1.inART.inc.comb <- dataprep.absnum(S1.inART.inc)
S1.inART.ir.comb <- dataprep.absnum(S1.inART.ir)
S1.inART.sir.comb <- dataprep.absnum(S1.inART.sir)
S1.inART.hiv.comb <- dataprep.hivprev(S1.inART.hiv)

S2.inART.inc.comb <- dataprep.absnum(S2.inART.inc)
S2.inART.ir.comb <- dataprep.absnum(S2.inART.ir)
S2.inART.sir.comb <- dataprep.absnum(S2.inART.sir)
S2.inART.hiv.comb <- dataprep.hivprev(S2.inART.hiv)

S3.inART.inc.comb <- dataprep.absnum(S3.inART.inc)
S3.inART.ir.comb <- dataprep.absnum(S3.inART.ir)
S3.inART.sir.comb <- dataprep.absnum(S3.inART.sir)
S3.inART.hiv.comb <- dataprep.hivprev(S3.inART.hiv)


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
S0b.stART.ir.comb$scenario <- "Baseline cytology and vaccination, no ART scale-up"
S0.inART.ir.comb$scenario <- "Baseline cytology and vaccination, with ART scale-up"
S1.inART.ir.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.inART.ir.comb$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
S3.inART.ir.comb$scenario <- "S2 + more frequent screening for WLHIV"

ir.comb.inART <- rbind.data.frame(S0b.stART.ir.comb, S0.inART.ir.comb, S1.inART.ir.comb, S2.inART.ir.comb, S3.inART.ir.comb)
ir.comb.inART$scenario <- factor(ir.comb.inART$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                        "Baseline cytology and vaccination, with ART scale-up",
                                                        "Scaled up HPV testing and 90% vaccination",
                                                        "S1 + 50% catch-up vaccination for WLHIV",
                                                        "S2 + more frequent screening for WLHIV"))
ir.comb.inART$outcome <- "Crude Incidence Rates"

# Combine into a DF for standardized IR
S0b.stART.sir.comb$scenario <- "Baseline cytology and vaccination, no ART scale-up"
S0.inART.sir.comb$scenario <- "Baseline cytology and vaccination, with ART scale-up"
S1.inART.sir.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.inART.sir.comb$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
S3.inART.sir.comb$scenario <- "S2 + more frequent screening for WLHIV"


sir.comb.inART <- rbind.data.frame(S0b.stART.sir.comb, S0.inART.sir.comb, S1.inART.sir.comb, S2.inART.sir.comb, S3.inART.sir.comb)
sir.comb.inART$scenario <- factor(sir.comb.inART$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                                    "Baseline cytology and vaccination, with ART scale-up",
                                                                    "Scaled up HPV testing and 90% vaccination",
                                                                    "S1 + 50% catch-up vaccination for WLHIV",
                                                                    "S2 + more frequent screening for WLHIV"))
sir.comb.inART$outcome <- "AS Incidence Rates"


##########
# Combine incident case counts for plotting # ----
##########

# Combine into a DF
S0b.stART.inc.comb$scenario <- "Baseline cytology and vaccination, no ART scale-up"
S0.inART.inc.comb$scenario <- "Baseline cytology and vaccination, with ART scale-up"
S1.inART.inc.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.inART.inc.comb$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
S3.inART.inc.comb$scenario <- "S2 + more frequent screening for WLHIV"

inc.comb.inART <- rbind.data.frame(S0b.stART.inc.comb, S0.inART.inc.comb, S1.inART.inc.comb, S2.inART.inc.comb, S3.inART.inc.comb)
inc.comb.inART$scenario <- factor(inc.comb.inART$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                        "Baseline cytology and vaccination, with ART scale-up",
                                                        "Scaled up HPV testing and 90% vaccination",
                                                        "S1 + 50% catch-up vaccination for WLHIV",
                                                        "S2 + more frequent screening for WLHIV"))
inc.comb.inART$outcome <- "Incident cases"


ratesandcases.inART <- rbind.data.frame(inc.comb.inART, ir.comb.inART, sir.comb.inART)

##########
# Combine HIV prevalence across scenarios # ----
##########

# Combine into a DF
S0b.stART.hiv.comb$scenario <- "Baseline cytology and vaccination, no ART scale-up"
S0.inART.hiv.comb$scenario <- "Baseline cytology and vaccination, with ART scale-up"
S1.inART.hiv.comb$scenario <- "Scaled up HPV testing and 90% vaccination"
S2.inART.hiv.comb$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
S3.inART.hiv.comb$scenario <- "S2 + more frequent screening for WLHIV"

hiv.comb.inART <- rbind.data.frame(S0b.stART.hiv.comb, S0.inART.hiv.comb, S1.inART.hiv.comb, S2.inART.hiv.comb, S3.inART.hiv.comb)
hiv.comb.inART$scenario <- factor(hiv.comb.inART$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                          "Baseline cytology and vaccination, with ART scale-up",
                                                          "Scaled up HPV testing and 90% vaccination",
                                                          "S1 + 50% catch-up vaccination for WLHIV",
                                                          "S2 + more frequent screening for WLHIV"))
hiv.comb.inART$outcome <- "HIV prevalence"



##########
# CALCULATE PAF in each year and cumulative PAF for each scenario # ----
##########

(prophiv.S0b.stART <- calc_prophiv(S0b.stART.inc.comb))
(ccHivRatio.S0b.stART <- calc_ccHivRatio(prophiv.S0b.stART , S0b.stART.hiv.comb))
(prophivVS.S0b.stART <- calc_prophivVS(S0b.stART.inc.comb))
(prophivVSofHIV.S0b.stART <- calc_prophivVSofHIV(S0b.stART.inc.comb))
calc_prophiv_cuml(S0b.stART.inc.comb, 2021, 2121)

(prophiv.S0.inART <- calc_prophiv(S0.inART.inc.comb))
(ccHivRatio.S0.inART <- calc_ccHivRatio(prophiv.S0.inART , S0.inART.hiv.comb))
(prophivVS.S0.inART <- calc_prophivVS(S0.inART.inc.comb))
(prophivVSofHIV.S0.inART <- calc_prophivVSofHIV(S0.inART.inc.comb))
calc_prophiv_cuml(S0.inART.inc.comb, 2021, 2121)

(prophiv.S1.inART <- calc_prophiv(S1.inART.inc.comb))
(ccHivRatio.S1.inART <- calc_ccHivRatio(prophiv.S1.inART , S1.inART.hiv.comb))
(prophivVS.S1.inART <- calc_prophivVS(S1.inART.inc.comb))
(prophivVSofHIV.S1.inART <- calc_prophivVSofHIV(S1.inART.inc.comb))
calc_prophiv_cuml(S1.inART.inc.comb, 2021, 2121)

(prophiv.S2.inART <- calc_prophiv(S2.inART.inc.comb))
(ccHivRatio.S2.inART <- calc_ccHivRatio(prophiv.S2.inART , S2.inART.hiv.comb))
(prophivVS.S2.inART <- calc_prophivVS(S2.inART.inc.comb))
(prophivVSofHIV.S2.inART <- calc_prophivVSofHIV(S2.inART.inc.comb))
calc_prophiv_cuml(S2.inART.inc.comb, 2021, 2121)

(prophiv.S3.inART <- calc_prophiv(S3.inART.inc.comb))
(ccHivRatio.S3.inART <- calc_ccHivRatio(prophiv.S3.inART , S3.inART.hiv.comb))
(prophivVS.S3.inART <- calc_prophivVS(S3.inART.inc.comb))
(prophivVSofHIV.S3.inART <- calc_prophivVSofHIV(S3.inART.inc.comb))
calc_prophiv_cuml(S3.inART.inc.comb, 2021, 2121)

# Combine into a data frame
prophiv.S0b.stART$scenario <- "Baseline cytology and vaccination, no ART scale-up"
prophiv.S0.inART$scenario <- "Baseline cytology and vaccination, with ART scale-up"
prophiv.S1.inART$scenario <- "Scaled up HPV testing and 90% vaccination"
prophiv.S2.inART$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
prophiv.S3.inART$scenario <- "S2 + more frequent screening for WLHIV"

prophiv.inART.comb <- rbind.data.frame(prophiv.S0b.stART, prophiv.S0.inART, prophiv.S1.inART, prophiv.S2.inART, prophiv.S3.inART)
prophiv.inART.comb$scenario <- factor(prophiv.inART.comb$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                          "Baseline cytology and vaccination, with ART scale-up",
                                                          "Scaled up HPV testing and 90% vaccination",
                                                          "S1 + 50% catch-up vaccination for WLHIV",
                                                          "S2 + more frequent screening for WLHIV"))

ccHivRatio.S0b.stART$scenario <- "Baseline cytology and vaccination, no ART scale-up"
ccHivRatio.S0.inART$scenario <- "Baseline cytology and vaccination, with ART scale-up"
ccHivRatio.S1.inART$scenario <- "Scaled up HPV testing and 90% vaccination"
ccHivRatio.S2.inART$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
ccHivRatio.S3.inART$scenario <- "S2 + more frequent screening for WLHIV"

ccHivRatio.inART.comb <- rbind.data.frame(ccHivRatio.S0b.stART, ccHivRatio.S0.inART, ccHivRatio.S1.inART, ccHivRatio.S2.inART, ccHivRatio.S3.inART)
ccHivRatio.inART.comb$scenario <- factor(ccHivRatio.inART.comb$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                                              "Baseline cytology and vaccination, with ART scale-up",
                                                                              "Scaled up HPV testing and 90% vaccination",
                                                                              "S1 + 50% catch-up vaccination for WLHIV",
                                                                              "S2 + more frequent screening for WLHIV"))

prophivVS.S0b.stART$scenario <- "Baseline cytology and vaccination, no ART scale-up"
prophivVS.S0.inART$scenario <- "Baseline cytology and vaccination, with ART scale-up"
prophivVS.S1.inART$scenario <- "Scaled up HPV testing and 90% vaccination"
prophivVS.S2.inART$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
prophivVS.S3.inART$scenario <- "S2 + more frequent screening for WLHIV"

prophivVS.inART.comb <- rbind.data.frame(prophivVS.S0b.stART, prophivVS.S0.inART, prophivVS.S1.inART, prophivVS.S2.inART, prophivVS.S3.inART)
prophivVS.inART.comb$scenario <- factor(prophivVS.inART.comb$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                                              "Baseline cytology and vaccination, with ART scale-up",
                                                                              "Scaled up HPV testing and 90% vaccination",
                                                                              "S1 + 50% catch-up vaccination for WLHIV",
                                                                              "S2 + more frequent screening for WLHIV"))

prophivVSofHIV.S0b.stART$scenario <- "Baseline cytology and vaccination, no ART scale-up"
prophivVSofHIV.S0.inART$scenario <- "Baseline cytology and vaccination, with ART scale-up"
prophivVSofHIV.S1.inART$scenario <- "Scaled up HPV testing and 90% vaccination"
prophivVSofHIV.S2.inART$scenario <- "S1 + 50% catch-up vaccination for WLHIV"
prophivVSofHIV.S3.inART$scenario <- "S2 + more frequent screening for WLHIV"

prophivVSofHIV.inART.comb <- rbind.data.frame(prophivVSofHIV.S0b.stART, prophivVSofHIV.S0.inART, prophivVSofHIV.S1.inART, prophivVSofHIV.S2.inART, prophivVSofHIV.S3.inART)
prophivVSofHIV.inART.comb$scenario <- factor(prophivVSofHIV.inART.comb$scenario, levels = c("Baseline cytology and vaccination, no ART scale-up", 
                                                                                  "Baseline cytology and vaccination, with ART scale-up",
                                                                                  "Scaled up HPV testing and 90% vaccination",
                                                                                  "S1 + 50% catch-up vaccination for WLHIV",
                                                                                  "S2 + more frequent screening for WLHIV"))



############################
## PLOTS
############################

colors <- RColorBrewer::brewer.pal(5, "Set2")

## Look at *median* outcomes across scenarios.
sir.comb.inART %>%  # Change the name of the dataframe to different outcomes (ir.comb.inART, sir.comb.inART, inc.comb.inART, )
  filter(group %in% c("Total", "HIVneg", "HIVpos_all" , "HIVpos_noART", "HIVpos_VS")) %>%
  ggplot() + 
    geom_line(aes(x=year, y = median, colour = scenario)) +
    facet_wrap(vars(group), scales = "free")

ggplot(prophiv.inART.comb) + 
    geom_line(aes(x=year, y = median, colour = scenario)) +
    #geom_line(data = hiv.comb.inART, aes(x=year, y = median, colour = scenario)) +
    geom_line(data = hiv.comb.inART, aes(x=year, y = median, colour = scenario), linetype = "dashed")


### Incidence rates
## Building up the scenarios for presentation

# Scenario 0b - Total pop only
png("sirS0b_tot.png", width = 800, height = 400)
#pdf("sirS0b_tot.pdf", width = 10, height = 5)
S0b.stART.sir.comb %>%
  filter(group == "Total") %>%
  ggplot() + 
    geom_line(aes(x = year, y = median), colour = "black", size = 1.25) +
    scale_y_continuous(breaks = seq(0, 105, 25), limits = c(0, 115)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Age-standardized incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black")) +
    plot_background 
dev.off()

# Scenario 0b - Adding in HIV states
png("sirS0b_byhiv.png", width = 800, height = 400)
#pdf("sirS0b_byhiv.pdf", width = 10, height = 5) # legend position for this size should be c(0.85, 0.85)
#pdf("sirS0b_byhiv.small.pdf", width = 6, height = 5)  # save it again at half size, and change legend position to c(0.8, 0.85)
S0b.stART.sir.comb %>%
  filter(group %in% c("Total", "HIVneg", "HIVpos_all")) %>%
  ggplot() + 
    geom_line(aes(x = year, y = median, colour = group), size = 1.25) +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos_all"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Age-standardized incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()

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

# Add in S0 with increasing ART
png("sirS0bS0_byhiv_inART.png", width = 800, height = 400)
#pdf("sirS0bS0_byhiv_inART.pdf", width = 6, height = 5)
ggplot() + 
    geom_line(data = filter(S0b.stART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.75, alpha = 0.45) +
    geom_line(data = filter(S0b.stART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1.5, alpha = 0.45) +
    geom_line(data = filter(S0b.stART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1.5, alpha = 0.45) +
    geom_line(data = filter(S0.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos_all"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
    #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
    labs(x = "Year", y = "Age-standardized incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()

png("irS0bS0_byhiv_inART.png", width = 800, height = 400)
#pdf("irS0bS0_byhiv_inART.pdf", width = 6, height = 5)
ggplot() + 
  geom_line(data = filter(S0b.stART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total", alpha = "artLevel"), size = 1.75) +
  geom_line(data = filter(S0b.stART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg", alpha = "artLevel"), size = 1.5) +
  geom_line(data = filter(S0b.stART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all", alpha = "artLevel"), size = 1.5) +
  geom_line(data = filter(S0.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total" , alpha = "artScale"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg" , alpha = "artScale"), size = 1) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all" , alpha = "artScale"), size = 1) +
  scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                      breaks = c("Total", "HIVneg", "HIVpos_all"),
                      labels = c("All women", "HIV-negative", "Women with HIV"),
                      name = "HIV status") +
  scale_alpha_manual(values = c("artLevel" = 0.45, "artScale" = 1),
                      breaks = c("artLevel", "artScale"),
                      labels = c("Baseline", "Baseline, with ART scale-up"),
                      name = "Scenario") +
  scale_y_continuous(breaks = seq(0, 250, 25), limits = c(0, 250)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "Year", y = "Crude cervical cancer incidence rate per 100,000 women") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.position = c(0.35, 0.85),
        legend.box = "horizontal") +
  plot_background 
dev.off()

# Add in S1 with increasing ART
#png("sirS0S1_byhiv.png", width = 800, height = 400)
pdf("sirS0S1_byhiv.pdf", width = 6, height = 5)
ggplot() + 
    geom_line(data = filter(S0.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
    geom_line(data = filter(S1.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dashed") +
    geom_line(data = filter(S1.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dashed") +
    geom_line(data = filter(S1.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dashed") +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos_all"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
    #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
    labs(x = "Year", y = "Age-standardized incidence rate per 100,000 women") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()

#png("irS0S1_byhiv.png", width = 800, height = 400)
pdf("irS0S1_byhiv.pdf", width = 6, height = 5)
ggplot() + 
  geom_line(data = filter(S0.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
  geom_line(data = filter(S1.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dashed") +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dashed") +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dashed") +
  scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                      breaks = c("Total", "HIVneg", "HIVpos_all"),
                      labels = c("All females", "HIV-negative", "HIV-positive"),
                      name = NULL) +
  scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "Year", y = "Age-standardized incidence rate per 100,000 women") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.8, 0.85)) +
  plot_background 
dev.off()


# Add in S2,S3 with increasing ART
#png("sirS0S1S2S3_byhiv.png", width = 800, height = 400)
#pdf("sirS0S1S2S3_byhiv.pdf", width = 10, height = 5)
pdf("sirS0S1S2S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
    geom_line(data = filter(S0.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
    geom_line(data = filter(S0.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
    geom_line(data = filter(S1.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dashed") +
    geom_line(data = filter(S1.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dashed") +
    geom_line(data = filter(S1.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dashed") +
    geom_line(data = filter(S2.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dotdash") +
    geom_line(data = filter(S2.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dotdash") +
    geom_line(data = filter(S2.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dotdash") +
    geom_line(data = filter(S3.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dotted") +
    geom_line(data = filter(S3.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dotted") +
    geom_line(data = filter(S3.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dotted") +
    scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                        breaks = c("Total", "HIVneg", "HIVpos_all"),
                        labels = c("All females", "HIV-negative", "HIV-positive"),
                        name = NULL) +
    scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "year", y = "Age-standardized incidence rate per 100,000 women") +
    theme(axis.title = element_text(size =12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.85)) +
    plot_background 
dev.off()

#png("irS0S1S2S3_byhiv.png", width = 800, height = 400)
#pdf("irS0S1S2S3_byhiv.pdf", width = 10, height = 5)
pdf("irS0S1S2S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
  geom_line(data = filter(S0.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
  geom_line(data = filter(S1.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dashed") +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dashed") +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dashed") +
  geom_line(data = filter(S2.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dotdash") +
  geom_line(data = filter(S2.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dotdash") +
  geom_line(data = filter(S2.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dotdash") +
  geom_line(data = filter(S3.inART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dotted") +
  geom_line(data = filter(S3.inART.ir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dotted") +
  geom_line(data = filter(S3.inART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dotted") +
  scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                      breaks = c("Total", "HIVneg", "HIVpos_all"),
                      labels = c("All females", "HIV-negative", "HIV-positive"),
                      name = NULL) +
  scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
  scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  labs(x = "year", y = "Crude incidence rate per 100,000 women") +
  theme(axis.title = element_text(size =12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.8, 0.85)) +
  plot_background 
dev.off()

# Add in S3 with increasing ART
#png("sirS0S1S3_byhiv.png", width = 800, height = 400)
#pdf("sirS0S1S3_byhiv.pdf", width = 10, height = 5)
pdf("sirS0S1S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
  geom_line(data = filter(S0.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25) +
  geom_line(data = filter(S0.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1) +
  geom_line(data = filter(S0.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1) +
  geom_line(data = filter(S1.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dashed") +
  geom_line(data = filter(S1.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dashed") +
  geom_line(data = filter(S1.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dashed") +
  geom_line(data = filter(S3.inART.sir.comb, group == "Total"), aes(x = year, y = median, colour = "Total"), size = 1.25, linetype = "dotted") +
  geom_line(data = filter(S3.inART.sir.comb, group == "HIVneg"), aes(x = year, y = median, colour = "HIVneg"), size = 1, linetype = "dotted") +
  geom_line(data = filter(S3.inART.sir.comb, group == "HIVpos_all"), aes(x = year, y = median, colour = "HIVpos_all"), size = 1, linetype = "dotted") +
  scale_colour_manual(values = c("Total" = "black", "HIVneg" = "#38b1b1", "HIVpos_all" = colors[3]),
                      breaks = c("Total", "HIVneg", "HIVpos_all"),
                      labels = c("All females", "HIV-negative", "HIV-positive"),
                      name = NULL) +
  scale_y_continuous(breaks = seq(0, 200, 25), limits = c(0, 200)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "year", y = "Age-standardized incidence rate per 100,000 women") +
  theme(axis.title = element_text(size =12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.8, 0.85)) +
  plot_background 
dev.off()

#################################################### FIGURE 1A #############################################################
# Crude ICC among all women
#png("irS0S1S3_byhiv.png", width = 800, height = 400)
#pdf("irS0S1S3_byhiv.pdf", width = 10, height = 5)
pdf("irS0S1S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
  geom_line(data = filter(S0b.stART.ir.comb, group == "Total"), aes(x = year, y = median, colour = "artLevel", linetype = "Total"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "Total", year > "2019"), aes(x = year, y = median, colour = "artScale", linetype = "Total"), size = 1.25) +
  geom_line(data = filter(S1.inART.ir.comb, group == "Total", year > "2019"), aes(x = year, y = median, colour = "hpvEnhanc", linetype = "Total"), size = 1.25) +
  geom_line(data = filter(S3.inART.ir.comb, group == "Total", year > "2019"), aes(x = year, y = median, colour = "hivEnhanc", linetype = "Total"), size = 1.25) +
  scale_colour_manual(values = c("artLevel" = "black", "artScale" = "grey45", "hpvEnhanc" = "dodgerblue2", "hivEnhanc" = "aquamarine3"),
                     breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                     labels = c("Baseline" , "ART scale-up only" , "Enhanced cervical cancer interventions" , "Enhanced cervical cancer interventions for women with HIV"),
                     name = "Scenario") +
  scale_linetype_manual(values = c("Total" = "solid", "HIVpos_all" = "dotted" , "HIVneg" = "dashed"),
                        breaks = c("Total", "HIVpos_all" , "HIVneg"),
                        labels = c("All women", "Women with HIV" , "Women without HIV"),
                        name = "HIV status") +
  scale_y_continuous(breaks = seq(0, 200, 50), limits = c(0, 200)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "", y = "Crude cervical cancer incidence rate per 100,000 women") +
  theme(axis.title = element_text(size =11),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1.5, 'cm'),
        legend.position = "bottom",
        legend.direction="vertical") +
  plot_background 
dev.off()

################################################## FIGURE 1b ##########################################################
# Crude ICC among HIV-negative and women with HIV
pdf("irS0S1S3_byhiv.small.pdf", width = 6, height = 5) # half size
ggplot() + 
  geom_line(data = filter(S0b.stART.ir.comb, group == "HIVpos_all"), aes(x = year, y = median, linetype = "HIVpos_all", colour = "artLevel"), size = 1.25) +
  geom_line(data = filter(S0b.stART.ir.comb, group == "HIVneg"), aes(x = year, y = median, linetype = "HIVneg", colour = "artLevel"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVpos_all", year > "2019"), aes(x = year, y = median, linetype = "HIVpos_all", colour = "artScale"), size = 1.25) +
  geom_line(data = filter(S0.inART.ir.comb, group == "HIVneg", year > "2019"), aes(x = year, y = median, linetype = "HIVneg", colour = "artScale"), size = 1.25) +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVpos_all", year > "2019"), aes(x = year, y = median, linetype = "HIVpos_all", colour = "hpvEnhanc"), size = 1.25) +
  geom_line(data = filter(S1.inART.ir.comb, group == "HIVneg", year > "2019"), aes(x = year, y = median, linetyper = "HIVneg", colour = "hpvEnhanc"), size = 1.25) +
  geom_line(data = filter(S3.inART.ir.comb, group == "HIVpos_all", year > "2019"), aes(x = year, y = median, linetype = "HIVpos_all", colour = "hivEnhanc"), size = 1.25) +
    geom_line(data = filter(S3.inART.ir.comb, group == "HIVneg", year > "2019"), aes(x = year, y = median, linetype = "HIVneg", colour = "hivEnhanc"), size = 1.25) +
  scale_linetype_manual(values = c("HIVpos_all" = "dotted" , "HIVneg" = "dashed"),
                      breaks = c("HIVpos_all" , "HIVneg"),
                      labels = c("Women with HIV" , "HIV-negative"),
                      name = "HIV status") +
  scale_colour_manual(values = c("artLevel" = "black", "artScale" = "grey45", "hpvEnhanc" = "dodgerblue2", "hivEnhanc" = "aquamarine3"),
                      breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                      labels = c("Baseline" , "ART scale-up only" , "Enhanced cervical cancer interventions" , "Enhanced cervical cancer interventions for women with HIV"),
                      name = "Scenario") +
  scale_y_continuous(breaks = seq(0, 200, 50), limits = c(0, 200)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "Year", y = "Crude cervical cancer incidence rate per 100,000 women") +
  theme(axis.title = element_text(size =11),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1.5, 'cm'),
        legend.position = "bottom",
        legend.direction="vertical") +
  plot_background 
dev.off()


### Proportion of cases in WLHIV and HIV prevalence

## Scenario S0b
pdf("prophiv.s0b.small.pdf", width = 6, height = 5) 
ggplot(prophiv.S0b.stART) +
    geom_area(aes(x = year, y = median), fill = colors[3]) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Proportion of cases among WLHIV") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

# Add HIV prevalence
pdf("prophiv.prev.s0b.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = prophiv.S0b.stART, aes(x = year, y = median), fill = colors[3]) +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") + 
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

## Scenario S0 with increasing ART
#pdf("prophiv.prev.s0bS0.inART.small.pdf", width = 6, height = 5) # save again at half size
#ggplot() +
#    geom_area(data = prophiv.S0b.stART, aes(x = year, y = median), fill = colors[3], alpha = 0.6) +
#    geom_area(data = prophiv.S0.inART, aes(x = year, y = median), fill = colors[3], alpha = 0.6) +
#    geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
#    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = median), fill = "#8B7DCA") + 
#    geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") + 
#    geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median), colour = "#3C2F74") +
#    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
#    #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
#    scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
#    labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
#    theme(axis.title = element_text(size = 12),
#          axis.text = element_text(size = 14),
#          axis.line = element_line(color = "black"),
#          legend.text = element_text(size = 14),
#          legend.position = c(0.85, 0.85)) +
#    plot_background
#dev.off()

## Scenario S0 with increasing ART
pdf("prophiv.prev.s0bS0.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophiv.S0b.stART, aes(x = year, y = median , colour = "propCases"), size = 1.75, alpha = 0.45) +
  geom_line(data = prophiv.S0.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25) +
  geom_line(data = S0b.stART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.75, alpha = 0.45) + 
  geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "Year", y = "Proportion") +
  scale_colour_manual(values = c("propCases" = "#486199", "hivPrev" = "#3C2F74"),
                      breaks = c("propCases", "hivPrev"),
                      labels = c("Proportion of cases in WLHIV", "HIV prevalence"),
                      name = NULL) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.85)) +
  plot_background
dev.off()

#Highlight the initial increase in PAF and then eventual decrease
S0b.min <- pmin(prophiv.S0b.stART$median, prophiv.S0.inART$median)
prophiv.S0b.min <- cbind.data.frame(year = prophiv.S0b.stART$year, min_prophiv = S0b.min)

pdf("prophiv.prev.s0bS0.inART.showincr.small.pdf", width = 6, height = 5) 
ggplot() +
    geom_area(data = prophiv.S0.inART, aes(x = year, y = median), fill = "#A4F283") +
    geom_area(data = prophiv.S0b.stART, aes(x = year, y = median), fill = colors[3], alpha = 0.6) +
    geom_area(data = prophiv.S0b.min, aes(x = year, y = min_prophiv), fill = colors[3]) +
    geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = median), fill = "#8B7DCA") + 
    geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") + 
    geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median), colour = "#3C2F74") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

pdf("prophiv.prev.s0bS0.inART.showdecr.small.pdf", width = 6, height = 5) 
ggplot() +
    geom_area(data = prophiv.S0.inART, aes(x = year, y = median), fill = "#A4F283") +
    geom_area(data = prophiv.S0b.stART, aes(x = year, y = median), fill = "#C76F98", alpha = 0.8) +
    geom_area(data = prophiv.S0b.min, aes(x = year, y = min_prophiv), fill = colors[3]) +
    geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
    geom_area(data = S0b.stART.hiv.comb, aes(x = year, y = median), fill = "#8B7DCA") + 
    geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") + 
    geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median), colour = "#3C2F74") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.85)) +
    plot_background
dev.off()

## Add in S1
pdf("prophiv.prev.s0s1.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
    geom_area(data = prophiv.S0.inART, aes(x = year, y = median, fill = "Baseline")) +
    geom_area(data = prophiv.S1.inART, aes(x = year, y = median, fill = "S1"), alpha = 0.55) +
    geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
    geom_line(data = prophiv.S1.inART, aes(x = year, y = median), colour = "#486199", linetype = "dashed") +
    # geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#738ABF") +
    # geom_line(data = prophiv.S1.inART, aes(x = year, y = median), colour = "#7AB4F5") +
    geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
    scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
    scale_fill_manual(values = c("Baseline" = colors[3], "S1" = "#97C4F7"),
                                          breaks = c("Baseline cytology and vaccination, with ART scale-up", "Scaled up HPV testing and 90% vaccination"),
                                          labels = c("Baseline cytology and vaccination, with ART scale-up", "Scaled up HPV testing and 90% vaccination"),
                                          name = NULL) +
    labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 14),
          legend.position = "none") +
    plot_background
dev.off()

## Add in S3
#pdf("prophiv.prev.s0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
#ggplot() +
#    geom_area(data = prophiv.S0.inART, aes(x = year, y = median, fill = "Baseline")) +
#    geom_area(data = prophiv.S1.inART, aes(x = year, y = median, fill = "S1"), alpha = 0.55) +
#    geom_area(data = prophiv.S3.inART, aes(x = year, y = median, fill = "S3")) +
#    geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
#    geom_line(data = prophiv.S1.inART, aes(x = year, y = median), colour = "#486199", linetype = "dashed") +
#    geom_line(data = prophiv.S3.inART, aes(x = year, y = median), colour = "#486199", linetype = "dotdash") +
#    geom_area(data = S1.inART.hiv.comb, aes(x = year, y = median), fill = "#A15AD1") + 
#    geom_area(data = S3.inART.hiv.comb, aes(x = year, y = median), fill = "#8771EB") +
#    geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") +
#    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
#    #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
#    scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
#    scale_fill_manual(values = c("Baseline" = colors[3], "S1" = "#97C4F7", "S3" = "#A0C9F8"),
#                      breaks = c("Baseline", "S1", "S3"),
#                      labels = c("Baseline cytology and vaccination, with ART scale-up", "Scaled up HPV testing and 90% vaccination",
#                                 "S2 + more frequent screening for WLHIV"),
#                      name = NULL) +
#    labs(x = "year", y = "Proportion of cases / HIV prevalence") +
#    theme(axis.title = element_text(size = 14),
#          axis.text = element_text(size = 14),
#          axis.line = element_line(color = "black"),
#          legend.text = element_text(size = 14),
#          legend.position = "none") +
#    plot_background
#dev.off()

## Add in S3
pdf("prophiv.prev.s0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophiv.S0.inART, aes(x = year, y = median, colour = "propCases"), size = 1.75) +
  geom_line(data = prophiv.S1.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25, linetype = "dashed") +
  geom_line(data = prophiv.S3.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25, linetype = "dotdash") +
  geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.75) + 
  geom_line(data = S1.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25, linetype = "dashed") +
  geom_line(data = S3.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25, linetype = "dotdash") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_colour_manual(values = c("propCases" = "#486199", "hivPrev" = "#3C2F74"),
                    breaks = c("propCases" , "hivPrev"),
                    labels = c("Proportion of cases in WLHIV", "HIV prevalence"),
                    name = NULL) +
  labs(x = "year", y = "Proportion of cases / HIV prevalence") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.85)) +
  plot_background
dev.off()

################################################# FIGURE 2 ###########################################################
## All scenarios - Proportion of cases and HIV prevalence
pdf("prophiv.prev.s0bs0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophiv.S0b.stART, aes(x = year, y = median, colour = "artLevel", linetype = "propCases"), size = 1.75) +
  geom_line(data = filter(prophiv.S0.inART, year > "2019"), aes(x = year, y = median, colour = "artScale", linetype = "propCases"), size = 1.75) +
  geom_line(data = filter(prophiv.S1.inART, year > "2019"), aes(x = year, y = median, colour = "hpvEnhanc", linetype = "propCases"), size = 1.75) +
  geom_line(data = filter(prophiv.S3.inART, year > "2019"), aes(x = year, y = median, colour = "hivEnhanc", linetype = "propCases"), size = 1.75) +
  geom_line(data = S0b.stART.hiv.comb, aes(x = year, y = median, colour = "artLevel", linetype = "hivPrev"), size = 1.75) + 
  geom_line(data = filter(S0.inART.hiv.comb, year > "2019"), aes(x = year, y = median, colour = "artScale", linetype = "hivPrev"), size = 1.75) + 
  geom_line(data = filter(S1.inART.hiv.comb, year > "2019"), aes(x = year, y = median, colour = "hpvEnhanc", linetype = "hivPrev"), size = 1.75) +
  geom_line(data = filter(S3.inART.hiv.comb, year > "2019"), aes(x = year, y = median, colour = "hivEnhanc", linetype = "hivPrev"), size = 1.75) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_colour_manual(values = c("artLevel" = "black", "artScale" = "grey45", "hpvEnhanc" = "dodgerblue2", "hivEnhanc" = "aquamarine3"),
                      breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                      labels = c("Baseline", "ART scale-up only" , "Enhanced cervical cancer interventions" , "Enhanced cervical cancer interventions for women with HIV"),
                      name = "Scenario") +
  scale_linetype_manual(values = c("propCases" = "solid", "hivPrev" = "dashed"),
                      breaks = c("propCases", "hivPrev"),
                      labels = c("Proportion cases among women with HIV" , "HIV prevalence among women"),
                      name = "Outcome") +
  labs(x = "Year", y = "Proportion") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1.5, 'cm'),
        legend.position = "bottom",
        legend.direction="vertical") +
  plot_background
dev.off()

## Add in S3 + S0b - HIV prevalence only
pdf("prophiv.prev.s0bs0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = S0b.stART.hiv.comb, aes(x = year, y = median, colour = "artLevel"), linetype = "dashed", size = 1.75, alpha = 0.45) + 
  geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median, colour = "artScale"), linetype = "dashed", size = 1.75) + 
  geom_line(data = S1.inART.hiv.comb, aes(x = year, y = median, colour = "hpvEnhanc"), linetype = "dashed", size = 1.25, linetype = "dashed") +
  geom_line(data = S3.inART.hiv.comb, aes(x = year, y = median, colour = "hivEnhanc"), linetype = "dashed", size = 1.25, linetype = "dotdash") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_colour_manual(values = c("artLevel" = "grey45", "artScale" = "black", "hpvEnhanc" = "dodgerblue2", "hivEnhanc" = "aquamarine3"),
                      breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                      labels = c("Baseline, no ART scale-up" , "Baseline" , "Enhanced HPV interventions" , "Enhanced HPV interventions for women with HIV"),
                      name = "Scenario") +
  labs(x = "Year", y = "HIV prevalence") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.position = c(0.35, 0.75)) +
  plot_background
dev.off()

## Add in S2
pdf("prophiv.prev.s0s1S2.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_area(data = prophiv.S0.inART, aes(x = year, y = median, fill = "Baseline")) +
  geom_area(data = prophiv.S1.inART, aes(x = year, y = median, fill = "S1"), alpha = 0.55) +
  geom_area(data = prophiv.S2.inART, aes(x = year, y = median, fill = "S3")) +
  geom_line(data = prophiv.S0.inART, aes(x = year, y = median), colour = "#486199") +
  geom_line(data = prophiv.S1.inART, aes(x = year, y = median), colour = "#486199", linetype = "dashed") +
  geom_line(data = prophiv.S2.inART, aes(x = year, y = median), colour = "#486199", linetype = "dotdash") +
  geom_area(data = S1.inART.hiv.comb, aes(x = year, y = median), fill = "#A15AD1") + 
  geom_area(data = S3.inART.hiv.comb, aes(x = year, y = median), fill = "#8771EB") +
  geom_area(data = S0.inART.hiv.comb, aes(x = year, y = median), fill = "#5F4BB6") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_fill_manual(values = c("Baseline" = colors[3], "S1" = "#97C4F7", "S3" = "#A0C9F8"),
                    breaks = c("Baseline", "S1", "S3"),
                    labels = c("Baseline cytology and vaccination, with ART scale-up", "Scaled up HPV testing and 90% vaccination",
                               "S1 + 50% catch-up vaccination for WLHIV"),
                    name = NULL) +
  labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = "none") +
  plot_background
dev.off()

# ggplot(paf.inART.comb) + 
#     geom_line(aes(x = year, y = paf_hiv, colour = scenario))
# ggplot(hiv.comb.inART) +
#     geom_line(aes(x = year, y = Prevalence, colour = scenario))







### Proportion of cases in WLHIV and WLHIV+VS and HIV prevalence

## Scenario S0 with increasing ART
pdf("prophivVS.prev.s0bS0.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophiv.S0b.stART, aes(x = year, y = median , colour = "propCases"), size = 1.75, alpha = 0.45) +
  geom_line(data = prophiv.S0.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25) +
  geom_line(data = prophivVS.S0b.stART, aes(x = year, y = median , colour = "propCasesVS"), size = 1.75, alpha = 0.45) +
  geom_line(data = prophivVS.S0.inART, aes(x = year, y = median, colour = "propCasesVS"), size = 1.25) +
  #geom_line(data = S0b.stART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.75, alpha = 0.45) + 
  #geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  labs(x = "Year", y = "Proportion") +
  scale_colour_manual(values = c("propCases" = "#486199", "propCasesVS" = "chartreuse3"), #, "hivPrev" = "#3C2F74"),
                      breaks = c("propCases", "propCasesVS"), # , "hivPrev"),
                      labels = c("Proportion of cases in WLHIV", "Proportion of cases in WLHIV+VS"), # , "HIV prevalence"),
                      name = NULL) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.85)) +
  plot_background
dev.off()

## Add in S3
pdf("prophivVS.prev.s0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophiv.S0.inART, aes(x = year, y = median, colour = "propCases"), size = 1.75) +
  geom_line(data = prophiv.S1.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25, linetype = "dashed") +
  geom_line(data = prophiv.S3.inART, aes(x = year, y = median, colour = "propCases"), size = 1.25, linetype = "dotdash") +
  geom_line(data = prophivVS.S0.inART, aes(x = year, y = median, colour = "propCasesVS"), size = 1.75) +
  geom_line(data = prophivVS.S1.inART, aes(x = year, y = median, colour = "propCasesVS"), size = 1.25, linetype = "dashed") +
  geom_line(data = prophivVS.S3.inART, aes(x = year, y = median, colour = "propCasesVS"), size = 1.25, linetype = "dotdash") +
  #geom_line(data = S1.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.75) + 
  #geom_line(data = S3.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25, linetype = "dashed") +
  #geom_line(data = S0.inART.hiv.comb, aes(x = year, y = median, colour = "hivPrev"), size = 1.25, linetype = "dotdash") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_colour_manual(values = c("propCases" = "#486199", "propCasesVS" = "chartreuse3"), #, "hivPrev" = "#3C2F74"),
                      breaks = c("propCases", "propCasesVS" ), #, "hivPrev"),
                      labels = c("Proportion of cases in WLHIV", "Proportion of cases in WLHIV+VS"), # , "HIV prevalence"),
                      name = NULL) +
  labs(x = "Year", y = "Proportion of cases / HIV prevalence") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.85)) +
  plot_background
dev.off()


### Ratio of proportion of cases within WLHIV to HIV prevalence

## S0, S0b, S1, S3
pdf("ccHivRatio.s0S0bs1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = ccHivRatio.S0b.stART, aes(x = year, y = median, linetype = "artLevel"), colour = "black", size = 1.75, alpha = 0.45) +
  geom_line(data = ccHivRatio.S0.inART, aes(x = year, y = median, linetype = "artScale"), colour = "black", size = 1.75) +
  geom_line(data = ccHivRatio.S1.inART, aes(x = year, y = median, linetype = "hpvEnhanc"), colour = "black", size = 1.25) +
  geom_line(data = ccHivRatio.S3.inART, aes(x = year, y = median, linetype = "hivEnhanc"), colour = "black", size = 1.25) +
  scale_y_continuous(breaks = seq(1, 5, 1), limits = c(1, 5)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_linetype_manual(values = c("artLevel" = "solid", "artScale" = "solid", "hpvEnhanc" = "dashed", "hivEnhanc" = "dotted"),
                        breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                        labels = c("Baseline" , "Baseline, with ART scale-up" , "Enhanced HPV interventions" , "Enhanced HPV interventions for women with HIV"),
                        name = "Scenario") +
  labs(x = "Year", y = "Ratio of proportion cases to HIV prevalence") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.position = c(0.35, 0.78)) +
  plot_background
dev.off()

ccHivRatio


################################################# FIGURE 4 ###########################################################
## All scenarios - Proportion of cases in WLHIV among WLHIV with VS and proportion VS
pdf("prophiv.prev.s0bs0s1S3.inART.small.pdf", width = 6, height = 5) # save again at half size
ggplot() +
  geom_line(data = prophivVSofHIV.S0b.stART, aes(x = year, y = median, colour = "artLevel"), size = 1.75) +
  geom_line(data = filter(prophivVSofHIV.S0.inART, year > "2019"), aes(x = year, y = median, colour = "artScale"), size = 1.75) +
  geom_line(data = filter(prophivVSofHIV.S1.inART, year > "2019"), aes(x = year, y = median, colour = "hpvEnhanc"), size = 1.75) +
  geom_line(data = filter(prophivVSofHIV.S3.inART, year > "2019"), aes(x = year, y = median, colour = "hivEnhanc"), size = 1.75) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  #scale_x_continuous(breaks = seq(2000 , 2120 , 20), limits = c(2000, 2120)) +
  scale_x_continuous(breaks = seq(2000 , 2070 , 10), limits = c(2000, 2071)) +
  scale_colour_manual(values = c("artLevel" = "black", "artScale" = "grey45", "hpvEnhanc" = "dodgerblue2", "hivEnhanc" = "aquamarine3"),
                      breaks = c("artLevel" , "artScale" , "hpvEnhanc" , "hivEnhanc"),
                      labels = c("Baseline", "ART scale-up only" , "Enhanced cervical cancer interventions" , "Enhanced cervical cancer interventions for women with HIV"),
                      name = "Scenario") +
  #scale_linetype_manual(values = c("propCases" = "solid"),
  #                      breaks = c("propCases"),
  #                      labels = c("Proportion of cases virally suppressed among women with HIV"),
  #                      name = "Outcome") +
  labs(x = "Year", y = "Proportion") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1.5, 'cm'),
        legend.position = "bottom",
        legend.direction="vertical") +
  plot_background
dev.off()



    
## Extract statistics
# Incidence rates in 2001, 2021, 2031, 2056, and 2071
#(startend.inART <- filter(S0.inART.ir.comb, year %in% c(2001, 2021, 2031, 2056, 2071), scenario == "Baseline cytology and vaccination, with ART scale-up"))
# Percent change in crude cervical cancer incidence rate over time by HIV status
#startend.inART %>% group_by(scenario) %>% mutate(pctchange_gen_2020 = Total[year == 2020] / Total[year == 2000],
#                                                 pctchange_gen_2070 = Total[year == 2070] / Total[year == 2020],
#                                                 pctchange_neg_2020 = HIVneg[year == 2020] / HIVneg[year == 2000],
#                                                 pctchange_neg_2070 = HIVneg[year == 2070] / HIVneg[year == 2020],
#                                                 pctchange_pos_2020 = HIVpos[year == 2020] / HIVpos[year == 2000],
#                                                 pctchange_pos_2070 = HIVpos[year == 2070] / HIVpos[year == 2020]) %>%
#    filter(year == 2000) %>%
#    select(scenario, pctchange_gen_2020, pctchange_gen_2070, pctchange_neg_2020, pctchange_neg_2070, pctchange_pos_2020, pctchange_pos_2070)
    
# Additional scenarios: percent change in crude cervical cancer incidence rate over time by HIV status
(startend.inART <- filter(ir.comb.inART, year %in% c(2001, 2021, 2031, 2056, 2071), scenario %in% c("Baseline cytology and vaccination, no ART scale-up",
                                                                                 "Baseline cytology and vaccination, with ART scale-up" ,                  
                                                                                 "Scaled up HPV testing and 90% vaccination",
                                                                                 "S2 + more frequent screening for WLHIV")))
#startend.inART %>% group_by(scenario) %>% mutate(pctchange_gen_2020 = Total[year == 2020] / Total[year == 2000],
#                                                 pctchange_gen_2070 = Total[year == 2070] / Total[year == 2020],
#                                           pctchange_neg_2020 = HIVneg[year == 2020] / HIVneg[year == 2000],
#                                           pctchange_neg_2070 = HIVneg[year == 2070] / HIVneg[year == 2020],
#                                           pctchange_pos_2020 = HIVpos[year == 2020] / HIVpos[year == 2000],
#                                           pctchange_pos_2070 = HIVpos[year == 2070] / HIVpos[year == 2020]) %>%
#    filter(year == 2000) %>%
#    select(scenario, pctchange_gen_2020, pctchange_gen_2070, pctchange_neg_2020, pctchange_neg_2070, pctchange_pos_2020, pctchange_pos_2070)
   
# Percent change in crude cervical cancer incidence relative to baseline scenario
pctreducS0S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducS0S0b) <- c('med' , 'min' , 'max')
pctreducS0S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducS0S0b[pctreducS0S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.inART[
       startend.inART$group == "Total" &
         startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
         startend.inART$year == i, c(5:29)] / 
         startend.inART[startend.inART$group == "Total" &
                          startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                          startend.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.inART[
         startend.inART$group == "Total" &
           startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
           startend.inART$year == i, c(5:29)] / 
           startend.inART[startend.inART$group == "Total" &
                            startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                            startend.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.inART[
         startend.inART$group == "Total" &
           startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
           startend.inART$year == i, c(5:29)] / 
           startend.inART[startend.inART$group == "Total" &
                            startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                            startend.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducS1S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducS1S0b) <- c('med' , 'min' , 'max')
pctreducS1S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducS1S0b[pctreducS1S0b$year == i , c(1:3)] <-
      c(apply((-1 + (startend.inART[
         startend.inART$group == "Total" &
         startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.inART$year == i, c(5:29)] / 
         startend.inART[startend.inART$group == "Total" &
         startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.inART$year == i, c(5:29)]))*100 , 1 , median) ,
      apply((-1 + (startend.inART[
         startend.inART$group == "Total" &
         startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.inART$year == i, c(5:29)] / 
         startend.inART[startend.inART$group == "Total" &
         startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.inART$year == i, c(5:29)]))*100 , 1 , min) , 
      apply((-1 + (startend.inART[
        startend.inART$group == "Total" &
        startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
        startend.inART$year == i, c(5:29)] / 
        startend.inART[startend.inART$group == "Total" &
        startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
        startend.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducS3S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducS3S0b) <- c('med' , 'min' , 'max')
pctreducS3S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducS3S0b[pctreducS3S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.inART[
       startend.inART$group == "Total" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "Total" &
       startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
       startend.inART$year == i, c(5:29)]))*100 , 1 , median) ,
     apply((-1 + (startend.inART[
       startend.inART$group == "Total" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "Total" &
       startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
       startend.inART$year == i, c(5:29)]))*100 , 1 , min) , 
     apply((-1 + (startend.inART[
       startend.inART$group == "Total" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "Total" &
       startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
       startend.inART$year == i, c(5:29)]))*100 , 1 , max))) }



# Additional scenarios: HIV prevalence over time
(startend.hivPrev.inART <- filter(hiv.comb.inART, year %in% c(2001, 2021, 2031, 2056, 2071), scenario %in% c("Baseline cytology and vaccination, no ART scale-up",
                                                                                                    "Baseline cytology and vaccination, with ART scale-up" ,                  
                                                                                                    "Scaled up HPV testing and 90% vaccination",
                                                                                                    "S2 + more frequent screening for WLHIV")))
# Percent change in HIV prevalence relative to baseline scenario
pctreducHivS0S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducHivS0S0b) <- c('med' , 'min' , 'max')
pctreducHivS0S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducHivS0S0b[pctreducHivS0S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducHivS1S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducHivS1S0b) <- c('med' , 'min' , 'max')
pctreducHivS1S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducHivS1S0b[pctreducHivS1S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.hivPrev.inART[
       startend.hivPrev.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducHivS3S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducHivS3S0b) <- c('med' , 'min' , 'max')
pctreducHivS3S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducHivS3S0b[pctreducHivS3S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.hivPrev.inART[
       startend.hivPrev.inART$scenario == "S2 + more frequent screening for WLHIV" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "S2 + more frequent screening for WLHIV" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.hivPrev.inART[
         startend.hivPrev.inART$scenario == "S2 + more frequent screening for WLHIV" & 
         startend.hivPrev.inART$year == i, c(5:29)] / 
         startend.hivPrev.inART[startend.hivPrev.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.hivPrev.inART$year == i, c(5:29)]))*100 , 1 , max))) }



# Additional scenarios: Proportion CC in women with HIV
(startend.prophiv.inART <- filter(prophiv.inART.comb, year %in% c(2001, 2021, 2031, 2056, 2071), scenario %in% c("Baseline cytology and vaccination, no ART scale-up",
                                                                                                             "Baseline cytology and vaccination, with ART scale-up" ,                  
                                                                                                             "Scaled up HPV testing and 90% vaccination",
                                                                                                             "S2 + more frequent screening for WLHIV")))
# Percent change in the proportion of CC in women with HIV relative to baseline scenario
pctreducPropHivS0S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducPropHivS0S0b) <- c('med' , 'min' , 'max')
pctreducPropHivS0S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducPropHivS0S0b[pctreducPropHivS0S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.prophiv.inART[
       startend.prophiv.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
         startend.prophiv.inART$year == i, c(5:29)] / 
         startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                  startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducPropHivS1S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducPropHivS1S0b) <- c('med' , 'min' , 'max')
pctreducPropHivS1S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducPropHivS1S0b[pctreducPropHivS1S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.prophiv.inART[
       startend.prophiv.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
         startend.prophiv.inART$year == i, c(5:29)] / 
         startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                  startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , max))) }

pctreducPropHivS3S0b <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(pctreducPropHivS3S0b) <- c('med' , 'min' , 'max')
pctreducPropHivS3S0b$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (pctreducPropHivS3S0b[pctreducPropHivS3S0b$year == i , c(1:3)] <-
     c(apply((-1 + (startend.prophiv.inART[
       startend.prophiv.inART$scenario == "S2 + more frequent screening for WLHIV" & 
         startend.prophiv.inART$year == i, c(5:29)] / 
         startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                  startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , median) ,
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "S2 + more frequent screening for WLHIV" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , min) , 
       apply((-1 + (startend.prophiv.inART[
         startend.prophiv.inART$scenario == "S2 + more frequent screening for WLHIV" & 
           startend.prophiv.inART$year == i, c(5:29)] / 
           startend.prophiv.inART[startend.prophiv.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                                    startend.prophiv.inART$year == i, c(5:29)]))*100 , 1 , max))) }



# Rate ratio of cervical cancer incidence among women with HIV compared to HIV-negative women
rateRatioS0b <- as.data.frame(matrix(ncol = 3, nrow = 5))
colnames(rateRatioS0b) <- c('med' , 'min' , 'max')
rateRatioS0b$year <- c(2001, 2021, 2031, 2056, 2071)
for(i in c(2001, 2021, 2031, 2056, 2071)) { 
  (rateRatioS0b[rateRatioS0b$year == i , c(1:3)] <-
     c(apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
         startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
         startend.inART$year == i, c(5:29)] / 
         startend.inART[startend.inART$group == "HIVneg" &
                          startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                          startend.inART$year == i, c(5:29)] , 1 , median) ,
       apply(startend.inART[
         startend.inART$group == "HIVpos_all" &
           startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
           startend.inART$year == i, c(5:29)] / 
           startend.inART[startend.inART$group == "HIVneg" &
                            startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                            startend.inART$year == i, c(5:29)] , 1 , min) , 
       apply(startend.inART[
         startend.inART$group == "HIVpos_all" &
           startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
           startend.inART$year == i, c(5:29)] / 
           startend.inART[startend.inART$group == "HIVneg" &
                            startend.inART$scenario == "Baseline cytology and vaccination, no ART scale-up" & 
                            startend.inART$year == i, c(5:29)] , 1 , max))) }

rateRatioS0 <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(rateRatioS0) <- c('med' , 'min' , 'max')
rateRatioS0$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (rateRatioS0[rateRatioS0$year == i , c(1:3)] <-
     c(apply((startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)]) , 1 , median) ,
     apply((startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)]) , 1 , min) , 
     apply((startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Baseline cytology and vaccination, with ART scale-up" & 
       startend.inART$year == i, c(5:29)]) , 1 , max))) }

rateRatioS1 <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(rateRatioS1) <- c('med' , 'min' , 'max')
rateRatioS1$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (rateRatioS1[rateRatioS1$year == i , c(1:3)] <-
     c(apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] , 1 , median) ,
     apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] , 1 , min) , 
     apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "Scaled up HPV testing and 90% vaccination" & 
       startend.inART$year == i, c(5:29)] , 1 , max))) }

rateRatioS3 <- as.data.frame(matrix(ncol = 3, nrow = 4))
colnames(rateRatioS3) <- c('med' , 'min' , 'max')
rateRatioS3$year <- c(2021, 2031, 2056, 2071)
for(i in c(2021, 2031, 2056, 2071)) { 
  (rateRatioS3[rateRatioS3$year == i , c(1:3)] <-
     c(apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] , 1 , median) ,
     apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] , 1 , min) , 
     apply(startend.inART[
       startend.inART$group == "HIVpos_all" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] / 
       startend.inART[startend.inART$group == "HIVneg" &
       startend.inART$scenario == "S2 + more frequent screening for WLHIV" & 
       startend.inART$year == i, c(5:29)] , 1 , max))) }

















#  !!!!!!!!!!!!!!!!!!!!!!NOT UPDATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        
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
