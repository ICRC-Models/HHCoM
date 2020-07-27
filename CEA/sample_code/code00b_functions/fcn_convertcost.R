# Subfunction conversions (same year, same country)
# scenario 0c: PPP given, LCU desired
# scenario 0d: USD given, LCU desired
# ... and reverse them:
# scenario 0a: LCU given, USD desired
# scenario 0b: LCU given, PPP desired

convert_lcu = function(val, iso3, yr, den, tofrom){
  tmp = wb(country = iso3, indicator = ifelse(den=="ppp", 'PA.NUS.PPP','PA.NUS.FCRF'), 
           mrv = 30, gapfill = TRUE)
  lcu_den = tmp$value[tmp$date==yr]
  
  return(val*ifelse(tofrom=="from", 1/lcu_den, lcu_den))
}

# change denomination (d2 != d2), same year
# scenario 0e: PPP given, USD desired
# scenario 0f: USD given, PPP desired

change_den = function(val, iso3, yr, d1, d2){
  if (any(c(d1, d2)=="lcu")){
    return(convert_lcu(val, iso3, yr, ifelse(d1=="lcu", d2, d1), ifelse(d1=="lcu", "from", "to")))
  } else {
    return(convert_lcu(convert_lcu(val, iso3, yr, d1, "to"), iso3, yr, d2, "from"))
  }
}

# Subfunction inflations
# scenario 1: LCU given, LCU desired
# others done by nesting
inflate_lcu = function(iso3, yr1, yr2){
  tmp = wb(country = iso3, indicator = 'FP.CPI.TOTL', mrv=30, gapfill = TRUE)
  return(tmp$value[tmp$date==yr2]/tmp$value[tmp$date==yr1])
}

# Big function
# Super-function, same country, different years same denominations.
# Super-function, same country, different years different denominations.
# Super-function, different countries, same year
# scenario 2a: USD given, USD desired, different countries
# scenario 2b: PPP given, PPP desired, different countries
# scenario 2c: LCU given, LCU desired, different countries
# nested... different countries, different years, potentially different denominations...

convertcost = function(val, t1, t2, c1, c2, d1, d2){
  
  # val: value given
  # t1: year of value cost
  # t2: year of value wanted
  # c1 is the country (ISO3) given as a string. 
  # d1 is the denomination given: "LCU", "PPP", "USD".
  # c2 is the country (ISO3) desired as a string. 
  # d2 is the denomination desired: "LCU", "PPP", "USD".
  
  # LATER: give option of providing the values of exchange, 
  # inflation, etc, but if blank, let it be WDI values
  
  # Values we will pull from WDI
  # PPP: 'PA.NUS.PPP' = "PPP conversion factor, GDP (LCU per international $)"
  # LCU per USD: 'PA.NUS.FCRF' = "Official exchange rate (LCU per US$, period average)"
  # GDP deflator in LCU: "NY.GDP.DEFL.KD.ZG.AD" = "Inflation, GDP deflator: linked series (annual %)"
  ## then use this: cumprod(rev(tmp$value)/100+1)*100
  # Use the GDP price deflator, as per WHO Guide to Cost-Effectiveness, section 3.2.6, page 43. 
  # Because of the volatility in the GDP price deflator, we are using the CPI
  
  # For Transferability of Prices across countries: 
  # See WHO Guide to Cost-Effectiveness, section 3.2.7, page 44
  # change to LCU in country 2 using PPPs in the year of the original data, then inflate.

if (t1==t2 & c1==c2){ # change denomination
  return(change_den(val, c1, t1, d1, d2))
} else if (c1==c2){ # change year & denomination if necessary
# convert to lcu, inflate, convert back in later year...
  tmp1 = ifelse(d1=="lcu", val, change_den(val, c1, t1, d1, "lcu"))
  tmp2 = tmp1*inflate_lcu(c1, t1, t2)
  return(ifelse(d2=="lcu", tmp2, change_den(tmp2, c1, t2, "lcu", d2)))
} else if (t1==t2){ # chance country and denomination (unless working in ppp)
  # convert to ppp (if not already), then convert to lcu of c2, and then to d2 (if not lcu)
  tmp1 = ifelse(d1=="ppp", val, change_den(val, c1, t1, d1, "ppp"))
  return(ifelse(d2=="ppp", tmp1, change_den(val, c2, t1, "ppp", d2)))
} else { # change country and year
  # convert to ppp (if not already), then convert to lcu of c2, then inflate, and then to d2
  tmp1 = ifelse(d1=="ppp", val, change_den(val, c1, t1, d1, "ppp"))
  tmp2 = change_den(tmp1, c2, t1, "ppp", "lcu")
  tmp3 = tmp1*inflate_lcu(c2, t1, t2)
  return(ifelse(d2=="lcu", tmp3, change_den(tmp2, c2, t2, "lcu", d2)))
}
 
}

