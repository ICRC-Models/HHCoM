# # Subfunction conversions (same year, same country)
# # scenario 0c: PPP given, LCU desired
# # scenario 0d: USD given, LCU desired
# # ... and reverse them:
# # scenario 0a: LCU given, USD desired
# # scenario 0b: LCU given, PPP desired
# 
# convert_lcu = function(val, iso3, yr, den, tofrom){
#   tmp = wb(country = iso3, indicator = ifelse(den=="ppp", 'PA.NUS.PPP','PA.NUS.FCRF'), mrv = 30, gapfill = TRUE)
#   lcu_den = tmp$value[tmp$date==yr]
#   
#   return(val*ifelse(tofrom=="from", lcu_den, 1/lcu_den))
# }
# # nested...
# # scenario 0e: PPP given, USD desired
# # scenario 0f: USD given, PPP desired
# 
# # Subfunctions inflations
# # scenario 1: LCU given, LCU desired
# # others: turn to LCU, then to the others...
# 
# inflate_lcu = function(val, iso3, yr1, yr2){
#   
#   tmp = wb(country = c2, indicator = 'FP.CPI.TOTL', 
#            startdate=min(yr1, yr2), enddate=max(yr1, yr2), gapfill = TRUE)
#   inf_yr1_yr2 = tmp$value[tmp$date==yr2]/tmp$value[tmp$date==yr1]
#   
# }
# 
# convertcost = function(val, t1, t2, c1, c2, d1, d2){
#   
#   # val: value given
#   # t1: year of value cost
#   # t2: year of value wanted
#   # c1 is the country (ISO3) given as a string. 
#   # d1 is the denomination given: "LCU", "PPP", "USD".
#   # c2 is the country (ISO3) desired as a string. 
#   # d2 is the denomination desired: "LCU", "PPP", "USD".
#   
#   # LATER: give option of providing the values of exchange, 
#   # inflation, etc, but if blank, let it be WDI values
#   
#   # Values we will pull from WDI
#   # PPP: 'PA.NUS.PPP' = "PPP conversion factor, GDP (LCU per international $)"
#   # LCU per USD: 'PA.NUS.FCRF' = "Official exchange rate (LCU per US$, period average)"
#   # GDP deflator in LCU: "NY.GDP.DEFL.KD.ZG.AD" = "Inflation, GDP deflator: linked series (annual %)"
#   ## then use this: cumprod(rev(tmp$value)/100+1)*100
#   # Use the GDP price deflator, as per WHO Guide to Cost-Effectiveness, section 3.2.6, page 43. 
#   # Because of the volatility in the GDP price deflator, we are using the CPI
#   
#   # For Transferability of Prices across countries: 
#   # See WHO Guide to Cost-Effectiveness, section 3.2.7, page 44
#   # change to LCU in country 2 using PPPs in the year of the original data, then inflate.
# 
# tmp = wb(country = c1, indicator = 'PA.NUS.PPP', mrv = 30, gapfill = TRUE)
# lcu_ppp_t1 = tmp$value[tmp$date==t1]
# lcu_ppp_t2 = tmp$value[tmp$date==t2]
# 
# tmp =  wb(country = c1, indicator = 'PA.NUS.FCRF', mrv = 30, gapfill = TRUE)
# lcu_usd_t1 = tmp$value[tmp$date==t1]
# lcu_usd_t2 = tmp$value[tmp$date==t2]
# 
# # For Transferability of Prices across countries: 
# # See WHO Guide to Cost-Effectiveness, section 3.2.7, page 44
# # change to LCU in country 2 using PPPs in the year of the original data, then inflate.
# 
# if (c1!=c2){
# tmp = wb(country = c2, indicator = 'PA.NUS.PPP', mrv = 10, gapfill = TRUE)
# lcu_ppp_t1_c2 = tmp$value[tmp$date==t1]
# 
# tmp =  wb(country = c2, indicator = 'PA.NUS.FCRF', mrv = 10, gapfill = TRUE)
# lcu_usd_t2_c2 = tmp$value[tmp$date==t2]
# }
# 
# # Inflation
# # tmp = wb(country = c1, indicator = "NY.GDP.DEFL.KD.ZG.AD", startdate=t1+1, enddate=t2, gapfill = F)
# # inf_drc_t1_t2 = cumprod(tmp1$value/100+1)
# tmp = wb(country = c2, indicator = 'FP.CPI.TOTL', startdate=t1+1, enddate=t2, gapfill = TRUE)
# inf_t1_t2 = tmp$value[tmp$date==t2]/tmp$value[tmp$date==t1]
# 
# 
# if ( & & ){
#   
# } else if{
#   
# }
# 
# 
# return(val*lcu_ppp_t1*inf_t1_t2/lcu_usd_t2)
# 
# }
# 
# # List of if...else statements
# 
# 
# 
# # Super-function, same country, different years same denominations.
# 
# # Super-function, same country, different years different denominations.
# 
# # Super-function, different countries, same year
# # scenario 2a: USD given, USD desired, different countries
# # scenario 2b: PPP given, PPP desired, different countries
# # scenario 2c: LCU given, LCU desired, different countries
# 
# # nested... different countries, different years, potentially different denominations...
# 
# # Print report in words... (as .rnw later)
# # including warnings when numbers had to be filled in with "gapfill"
# 
# 
# 
# 
# 
# 
