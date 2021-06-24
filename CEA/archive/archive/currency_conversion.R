###########################################
# Author: Miranda Tao, Bianca Zlavog
# Updated: 11/08/2019
#' Taking a formatted dataset and convert currency in values using IHME's deflator/xrates/ppp series
#' 
#' @param data [data.frame/data.table] A dataframe or datatable. MUST SPECIFY
#' @param col.loc [str] Column name that identify iso3 code, usually `ihme_loc_id`. MUST SPECIFY
#' @param col.value [str, vector] Column names that needs to be converted, take single or multiple. MUST SPECIFY
#' @param col.currency [str] Column that has the raw currency. This column takes values ['lcu','usd','eur','ppp'] and is not cap sensitive. LEAVE BLANK IF \param{currency} is specified
#' @param currency [str] Raw currency. LEAVE BLANK IF \param{col.currency} is specified
#' @param col.currency.year [str] Column that has raw currency year. LEAVE BLANK IF \param{currency.year} is specified
#' @param currency.year [int] Year of raw currency . LEAVE BLANK IF \param{col.currency.year} is specified
#' @param base.year  [int] Default is set to `the latest`. Year you want to base the output unit, default is 2017
#' @param base.unit [str] Default is `usd`. Please choose 'ppp' or 'usd'.
#' @param converter.version [num] Default is set to `the latest`. IHME Versions for xrates, deflator, PPP, choose from [1:2017/11, 2: 2018/5, 3: 2018/11, 3.1: 2019/2 fix VEN (final FGH 2018), 4.1: 2019/9 (initial estimates FGH 2019), 4.2: 2019/10 (updated IMF raw data), 4.3: 2019/11 (stop using WB XR data and old WHO PPP data), 4.4: 2019/11, (fix VEN), 4.5: 2019/11 (fix SOM and LBR XR), 4.6: 2019/12 (first submission FGH2019, fix PPP for VEN, STP, MRT, LBR), 4.7 (new WHO GHED data)]
#' @param inflate.usd [logical] Default is `FALSE`. If TRUE, nominal LCUs are converted to nominal USD, then inflated in USD
#' @param simplify  [logical] Default is `TRUE`. Return the simplified data with values converted, no other helper variables

#' @return [data.table] A datatable in the same format as input, but \param{col.value} converted





###########################################



### runtime configuration
if (Sys.info()["sysname"] %in% c("Linux","Darwin")){
  j <- "/home/j/" 
} else { 
  j <- "J:"
  if (!("pacman" %in% rownames(installed.packages()))){
    install.packages('pacman')
  }
}

require(pacman)
pacman::p_load(data.table, dplyr, tidyr, readstata13)


currency_conversion <- function(data,
                                col.loc, 
                                col.value, 
                                col.currency, 
                                currency,
                                col.currency.year, 
                                currency.year, 
                                base.year = 2019, 
                                base.unit = 'usd', 
                                converter.version = 4.7, 
                                inflate.usd = F, 
                                simplify = T 
                                
){
  
  
  data <- as.data.table(data)
  
  if (missing(col.currency)){
    col.currency <- 'currency'
    data[, eval(col.currency) := toupper(eval(currency))]
  }
  
  
  if (missing(col.currency.year)){
    col.currency.year <- 'currency_year'
    data[, eval(col.currency.year) := eval(currency.year)]
  }
  
  # force value types
  data[, eval(col.loc) := as.character(get(col.loc))]
  data[, eval(col.value) := lapply(eval(col.value), function(x) as.numeric(get(x)))]
  data[, eval(col.currency) := toupper(as.character(get(col.currency)))]
  data[, eval(col.currency.year) := as.integer(get(col.currency.year))]
  
  
  # create new names
  col.value.new = paste(col.value, '_new',sep = '')
  col.currency.new = paste0(col.currency, "_new")
  col.currency.year.new = paste0(col.currency.year, "_new")
  
  data[, eval(col.value.new) := lapply(eval(col.value), function(x) get(x))] # initialize converted value
  data[, eval(col.currency.new) := get(col.currency)] # initialize new currency
  data[,  eval(col.currency.year.new) := get(col.currency.year)] # initialize new currency year
  
  # create a list of columns to keep at the end for simplified version
  cols_keep <- names(data)
  cols_drop <- c(col.value, col.currency.new, col.currency.year.new,col.currency, col.currency.year)
  cols_keep <- cols_keep[!(cols_keep %in% cols_drop)]
  
  # check if there are NAs in location, value, currency or currency year
  if (nrow(data[get(col.loc) == "" | is.na(get(col.loc))]) != 0){
    warning(call. = F, 'There are NAs in location column')
  }
  
  # check if there are NAs in value
  if (nrow(data[is.na(get(col.value))]) != 0){
    warning(call. = F, 'There are NAs in value column')
  }
  
  # check if there are NAs in currency
  if (nrow(data[get(col.currency) == "" | is.na(get(col.currency))]) != 0){
    warning(call. = F, 'There are NAs in currency column')
  }
  
  # check if there are NAs in currency year
  if (nrow(data[is.na(get(col.currency.year))]) != 0){
    warning(call. = F, 'There are NAs in currency year column')
  }
  
  # check if there are non-standard currency
  if (nrow(data[!grepl('usd|eur|lcu|ppp',get(col.currency), ignore.case = T)]) != 0){
    warning(call. = F,paste0('Please check these non-standard currency eg: \"',
                             unique(data[!grepl('usd|eur|lcu',get(col.currency), ignore.case = T), get(col.currency)]), '\"'))
  }
  
  
  
  #============== get in all the dependent dataset =======================
  ## get eur exchange rate
  oecd_xrate = fread(paste0(j, "/Project/IRH/LDI_PPP/XR/data/OECD_XRATES_NattoUSD_1950_2019.csv"))
  oecd_xrate = oecd_xrate[LOCATION == "EA19", .(TIME, Value)]
  setnames(oecd_xrate, c('TIME','Value'), c('year_id', 'eur_usd'))
  oecd_xrate[, year_id:= as.integer(year_id)]
  
  #============== converter versions =======================
  if (converter.version == 1){
    
    x_rate_name = 'xr_output_20180516.dta'
    def_name = 'def_output_20171105.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/ppp_output_20171115.dta")) %>%  #newer: ppp_output_20180516.dta
      select(iso3, year,IHME_PPP_series) %>% data.table
    
    # PRK ppp
    prk_ppp <- read.dta13(paste0(j,"/Project/IRH/LDI_PPP/PPP/data/ppp_output_20180516.dta")) %>%
      select(iso3, year, IHME_PPP_series) %>% 
      filter(iso3 == 'PRK') %>% data.table
    ppp <- rbind(ppp, prk_ppp)
    
  } else if (converter.version == 2){
    
    x_rate_name = 'xr_output_20180516.dta'
    def_name = 'def_output_20180509.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/ppp_output_20180516.dta")) %>%  #newer: ppp_output_20180516.dta
      select(iso3, year,IHME_PPP_series) %>% data.table
    
  } else if (converter.version == 3){
    
    x_rate_name = 'xr_output_20180831 - VEN fix.dta'
    def_name = 'def_output_20180830.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/ppp_output_20181115 - VEN fix.dta")) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 3.1){
    
    x_rate_name = 'xr_output_20180831 - VEN fix.dta'
    def_name = 'def_output_20180830 - VEN fix 2.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/ppp_output_20181115 - VEN fix.dta")) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.1){
    
    x_rate_name = 'xr_output_20190930 - VEN fix.dta'
    def_name = 'def_output_20190923.dta'
    ppp_name = 'ppp_output_20190919 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.2){
    
    x_rate_name = 'xr_output_20191023 - VEN fix.dta'
    def_name = 'def_output_20191023.dta'
    ppp_name = 'ppp_output_20191023 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.3){
    
    x_rate_name = 'xr_output_20191029 - VEN fix.dta'
    def_name = 'def_output_20191023.dta'
    ppp_name = 'ppp_output_20191029 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.4){
    
    x_rate_name = 'xr_output_20191108 - VEN fix.dta'
    def_name = 'def_output_20191108 - VEN fix.dta'
    ppp_name = 'ppp_output_20191108 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.5){
    
    x_rate_name = 'xr_output_20191118 - VEN fix.dta'
    def_name = 'def_output_20191108 - VEN fix.dta'
    ppp_name = 'ppp_output_20191108 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.6){
    
    x_rate_name = 'xr_output_20191205 - VEN fix.dta'
    def_name = 'def_output_20191108 - VEN fix.dta'
    ppp_name = 'ppp_output_20191205 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  } else if (converter.version == 4.7){
    
    x_rate_name = 'xr_output_20191220 - VEN fix.dta'
    def_name = 'def_output_20191220 - VEN fix.dta'
    ppp_name = 'ppp_output_20191220 - VEN fix.dta'
    
    # get ppp series
    ppp = read.dta13(paste0(j, "/Project/IRH/LDI_PPP/PPP/data/",ppp_name)) %>%   
      dplyr::select(iso3, year,IHME_PPP_series) %>% data.table
  }
  
  ## get other exchange rates
  x_rate = as.data.table(read.dta13(paste0(j, "/Project/IRH/LDI_PPP/XR/data/",x_rate_name), generate.factors=T))
  setnames(x_rate, c('iso3','year','IHME_XR_series'), c('ihme_loc_id','year_id', 'lcu_usd'))
  x_rate = x_rate[,.(ihme_loc_id, year_id, lcu_usd)]
  
  # Fix formatting of year variable in this dataset
  if (converter.version == 4.1 | converter.version == 4.2 | converter.version == 4.3 | converter.version == 4.4 | converter.version == 4.5 | converter.version == 4.6 | converter.version == 4.7){
    x_rate <- x_rate[, year_id := as.numeric(as.character(year_id))]
  }  
  
  # get deflator series
  deflator = as.data.table(read.dta13(paste0(j, "/Project/IRH/LDI_PPP/DEF/data/",def_name)))
  setnames(deflator, c('iso3','year','IHME_DEF_series'), c('ihme_loc_id','year_id', 'deflator'))
  deflator = deflator[,.(ihme_loc_id, year_id, deflator)]
  
  
  
  #******************************************************
  ## need to rebase to x year you specify
  deflator_year = deflator[year_id == base.year, .(ihme_loc_id, deflator)]
  setnames(deflator_year, 'deflator', 'deflator_year')
  # if(unique(data$currency_year) != get(base.year)) {
  #   print(paste0("Rebasing deflators to ", base.year))
  # }
  deflator = deflator %>% 
    inner_join(deflator_year, by = 'ihme_loc_id') %>% 
    mutate(deflator = deflator/deflator_year) %>% 
    dplyr::select(-deflator_year) %>% data.table 
  
  # rename PPP series to be able to merge with the input data
  setnames(ppp, c('iso3', 'IHME_PPP_series', 'year'), c('ihme_loc_id','ppp','year_id'))
  
  #==============================================================================
  
  # split data based on input currency measure (ok if empty)
  eur = data[grepl('eur', get(col.currency), ignore.case = T), with = T]
  usd = data[grepl('usd', get(col.currency), ignore.case = T), with = T]
  lcu = data[grepl('lcu', get(col.currency), ignore.case = T), with = T]
  
  # if the data is already in base year PPP
  dt.ppp = data[grepl('ppp', get(col.currency), ignore.case = T), with = T]
  
  #=======================================
  # convert nominal euro to nominal usd
  #=======================================
  if (nrow(eur) != 0){
    eur <- merge(eur, oecd_xrate, 
                 by.x = col.currency.year.new, 
                 by.y = 'year_id', 
                 all.x = T)
    eur[, eval(col.value.new) := lapply(eval(col.value), function(x) get(x) / eur_usd)]
    eur[, eval(col.currency.new) := 'USD']
  }
  
  # combine nominal usd from usd subset and converted eur subset 
  usd <- rbindlist(list(usd, eur), fill = T, use.names = T)
  
  
  
  #############################################################################################################################################
  #                                                   IF OUTPUT IS IN BASE YEAR PPP
  #############################################################################################################################################
  
  if (base.unit %in% c('ppp','PPP')){
    
    
    #=======================================================
    # If we ONLY want to deflate between different PPP years
    #=======================================================
    if ((nrow(dt.ppp) != 0) & (unique(data[,get(col.currency)]) == 'PPP')){
      
      # step0: convert ppp to lcu
      dt.ppp <- merge(dt.ppp, ppp,
                      by.x=c(col.loc, col.currency.year.new),
                      by.y=c('ihme_loc_id', 'year_id'),
                      all.x = T)
      #print(paste0("Converting ", get(col.currency.year), " PPP to ", get(col.currency.year), " LCU"))
      dt.ppp[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)*ppp)]
      dt.ppp[, eval(col.currency.new)] = 'LCU'
      
      lcu <- rbindlist(list(lcu, dt.ppp), use.names = T, fill = T)
      
      # step2: inflate lcu to base year lcu
      lcu <- merge(lcu, deflator,
                   by.x=c(col.loc, col.currency.year.new),
                   by.y=c('ihme_loc_id', 'year_id'),
                   all.x = T)
      #print(paste0("Deflating ", get(col.currency.year), " LCU to ", get(base.year), " LCU"))
      lcu[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/deflator)]
      lcu[, eval(col.currency.year.new) := eval(base.year)]
      
      # step3: convert base year lcu to base year PPP
      data.converted <- merge(lcu, ppp,
                              by.x=c(col.loc, col.currency.year.new),
                              by.y=c('ihme_loc_id', 'year_id'),
                              all.x = T)
      #print(paste0("Converting ", get(base.year), " LCU to ", get(base.year), " PPP"))
      data.converted[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/ppp.y)]
      data.converted[, eval(col.currency.new)] = toupper(base.unit)
      
    } else if (nrow(dt.ppp) == 0){
      
      #-------------------------
      # if inflate.usd == FALSE
      #------------------------
      if (inflate.usd == F){
        
        #===========================================
        # if there is data in usd, convert to lcu
        #===========================================
        if (nrow(usd) != 0){
          # convert base year usd to base year lcu
          usd <- merge(usd, x_rate, 
                       by.x=c(col.loc, col.currency.year.new),
                       by.y=c('ihme_loc_id', 'year_id'), 
                       all.x = T)
          #print(paste0("Converting ", get(col.currency.year), " USD to ", get(col.currency.year), " LCU"))
          usd[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x) * lcu_usd)]
          usd[, eval(col.currency.new)] = 'LCU'
        }
        
        #================================
        # bind everything in LCU 
        #================================
        lcu <- rbindlist(list(lcu, usd), use.names = T, fill = T)
        
        # step1: inflate lcu to base year lcu
        lcu <- merge(lcu, deflator,
                     by.x=c(col.loc, col.currency.year.new),
                     by.y=c('ihme_loc_id', 'year_id'),
                     all.x = T)
        #print(paste0("Deflating ", get(col.currency.year), " LCU to ", get(base.year), " LCU"))
        lcu[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/deflator)]
        lcu[, eval(col.currency.year.new) := eval(base.year)]
        
        # step2: convert base year lcu dollar to base year PPP
        data.converted <- merge(lcu, ppp,
                                by.x=c(col.loc, col.currency.year.new),
                                by.y=c('ihme_loc_id', 'year_id'),
                                all.x = T)
        #print(paste0("Converting ", get(base.year), " LCU to ", get(base.year), " PPP"))
        data.converted[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/ppp)]
        data.converted[, eval(col.currency.new)] = toupper(base.unit)
        
      } else{
        #-------------------------------------------------
        # if inflate.usd == TRUE (WE DON'T DO IT THIS WAY)
        #-------------------------------------------------
        stop("You can't inflate in USD when converting from LCU/USD/EUR to PPP. Everything should be inflated in LCU.") 
        
      }
      
      
    } else{
      
      # !STOP! for any other conditions
      stop("Since you are starting with PPP, the dataset doesn't allow 'LCU' or 'USD' to coexist in the currency column. ")
    }
    
    
    #############################################################################################################################################
    #                                                 IF OUTPUT IS IN BASE YEAR USD
    #############################################################################################################################################
  } else if (base.unit %in% c('usd','USD')) { 
    
    #====================================================================================================================
    #                                           If all data are in x year real PPP
    #===================================================================================================================
    
    
    
    if (nrow(dt.ppp) != 0){
      if (unique(data[,get(col.currency)]) != 'PPP'){
        stop("Since you are starting with PPP, the dataset doesn't allow 'LCU' or 'USD' to coexist in the currency column. ")
      }
      
      # step0: convert ppp to lcu
      dt.ppp <- merge(dt.ppp, ppp,
                      by.x=c(col.loc, col.currency.year.new),
                      by.y=c('ihme_loc_id', 'year_id'),
                      all.x = T)
      #print(paste0("Converting ", get(col.currency.year), " PPP to ", get(col.currency.year), " LCU"))
      dt.ppp[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)*ppp)]
      dt.ppp[, eval(col.currency.new)] = 'LCU'
      
      lcu <- rbindlist(list(lcu, dt.ppp), use.names = T, fill = T)
      
      # step1: inflate lcu to base year lcu
      lcu <- merge(lcu, deflator,
                   by.x=c(col.loc, col.currency.year.new),
                   by.y=c('ihme_loc_id', 'year_id'),
                   all.x = T)
      #print(paste0("Deflating ", get(col.currency.year), " LCU to ", get(base.year), " LCU"))
      lcu[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/deflator)]
      lcu[, eval(col.currency.year.new) := eval(base.year)]
      
      
      # step 2: convert real lcu to real usd
      usd <- merge(lcu, x_rate,
                   by.x=c(col.loc, col.currency.year.new),
                   by.y=c('ihme_loc_id', 'year_id'),
                   all.x = T)
      #print(paste0("Converting ", get(base.year), " LCU to ", get(base.year), " USD"))
      usd[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x) / lcu_usd)]
      usd[, eval(col.currency.new)] = 'USD'
      
      data.converted <- copy(usd)
      
    } else{
      #==================================================================================================================
      #                                                If there is no data in PPP but in LCU/USD
      #==================================================================================================================
      
      #-------------------------
      # if inflate.usd == True
      #------------------------
      if (inflate.usd == T){
        
        #=====================================
        # exchange nominal LCU to nominal USD
        #=====================================
        if (nrow(lcu)!=0){
          lcu <- merge(lcu, x_rate, 
                       by.x=c(col.loc, col.currency.year.new),
                       by.y=c('ihme_loc_id', 'year_id'), 
                       all.x = T)
          #print(paste0("Converting ", get(col.currency.year), " LCU to ", get(col.currency.year), " USD"))
          lcu[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x) / lcu_usd)]
          lcu[, eval(col.currency.new)] = 'USD'
        }
        
        # bind everything in nominal USD
        usd <- rbindlist(list(lcu, usd), use.names = T, fill = T)
        
        
        #==================================
        # inflate nominal USD to real USD
        #==================================
        
        # get USA deflators
        deflator_usa <- deflator[ihme_loc_id == 'USA'][,ihme_loc_id:=NULL]
        usd <- merge(usd, deflator_usa,
                     by.x=c(col.currency.year.new),
                     by.y=c('year_id'),
                     all.x = T)
        #print(paste0("Deflating ", get(col.currency.year), " USD to ", get(base.year), " USD using USA deflators"))
        usd[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/deflator)]
        usd[, eval(col.currency.year.new) := eval(base.year)]
        
        data.converted <- copy(usd)
        
        
        
        #-------------------------
        # if inflate.usd == False
        #------------------------
      }else{
        #======================================================
        # exchange nominal USD to nominal LCU if any
        #======================================================
        if (nrow(usd)!=0){
          usd <- merge(usd, x_rate, 
                       by.x=c(col.loc, col.currency.year.new),
                       by.y=c('ihme_loc_id', 'year_id'), 
                       all.x = T)
          #print(paste0("Converting ", get(col.currency.year), " USD to ", get(col.currency.year), " LCU"))
          usd[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x) * lcu_usd)]
          usd[, eval(col.currency.new)] = 'LCU'
        }
        
        #================================================================
        # bind nominal LCU and converted nominal LCU, deflate to real LCU 
        #================================================================
        lcu <- rbindlist(list(lcu, usd), use.names = T, fill = T)
        
        # deflate nominal LCU to real LCU
        lcu <- merge(lcu, deflator,
                     by.x=c(col.loc, col.currency.year.new),
                     by.y=c('ihme_loc_id', 'year_id'),
                     all.x = T)
        #print(paste0("Deflating ", get(col.currency.year), " LCU to ", get(base.year), " LCU"))
        lcu[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x)/deflator)]
        lcu[, eval(col.currency.year.new) := eval(base.year)]
        
        
        #===============================================================
        # exchange real LCU to real USD
        #===============================================================
        if (nrow(lcu) != 0){
          
          # if xrate not being merged
          if ('lcu_usd' %in% names(lcu)) {
            lcu[,lcu_usd := NULL]
          }     
          
          # exchange real LCU to real USD
          lcu_to_usd <- merge(lcu, x_rate,
                              by.x=c(col.loc, col.currency.year.new),
                              by.y=c('ihme_loc_id', 'year_id'),
                              all.x = T)
          #print(paste0("Converting ", get(base.year), " LCU to ", get(base.year), " USD"))
          lcu_to_usd[, eval(col.value.new) := lapply(eval(col.value.new), function(x) get(x) / lcu_usd)]
          lcu_to_usd[, eval(col.currency.new)] = 'USD'
          
        }
        
        data.converted <- copy(lcu_to_usd)
        # data.converted[, eval(col.currency.new)] = toupper(base.unit)
        # data.converted[, eval(col.currency.year.new) := eval(base.year)]
      }
      
    }
    
    
  } else {
    stop(call. = F, 'the base unit you entered is not valid for this function, please choose from "usd", "ppp" ')
  }
  
  
  # check if there are NAs in value 
  if (NA %in% lcu$value_raw_new){
    warning(call. = F, 'There are NA values generated during conversions, please check')
  } 
  
  # simplify the output data
  if (simplify){
    data.converted <- data.converted[,c(cols_keep), with = F]
    setnames(data.converted, col.value.new, col.value)
  }
  
  return(data.converted)
}
