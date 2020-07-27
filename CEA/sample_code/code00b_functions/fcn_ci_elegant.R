#' Various functions to write summaries elegantly
ci_string = function(myvector){
  result = paste(round(sort(myvector)[2]), 
                 " (", round(min(myvector)), ", ", 
                 round(max(myvector)), ")", sep="")
  return(result)
}

mean_se_string = function(myvector){
  result = paste("Mean: ", round(sort(myvector)[2])/100, 
                 ", SE: ", round(max(myvector)-min(myvector))/4/100, sep="")
  return(result)
}

ci_string_1000 = function(myvector){
  if(min(myvector)>1000){
    result = paste(format(round(sort(myvector)[2]/1000, 0), nsmall=0, big.mark=","), " (", 
                   format(round(min(myvector)/1000, 0), nsmall=0, big.mark=","), ", ", 
                   format(round(max(myvector)/1000, 0), nsmall=0, big.mark=","), ")", sep="")
  }else{
    result = paste(format(round(sort(myvector)[2]/1000, 3), nsmall=0, big.mark=","), " (", 
                   format(round(min(myvector)/1000, 3), nsmall=0, big.mark=","), ", ", 
                   format(round(max(myvector)/1000, 3), nsmall=0, big.mark=","), ")", sep="")
  }
  return(result)
}

ci_string_comma = function(myvector){
  result = paste(format(round(sort(myvector)[2], 0), nsmall=0, big.mark=","), " (", 
                 format(round(min(myvector), 0), nsmall=0, big.mark=","), ", ", 
                 format(round(max(myvector), 0), nsmall=0, big.mark=","), ")", sep="")
  return(result)
}

options(scipen = 4)
options(digits = 4)


ci_string_dec = function(myvector, dec){
  # sign0 = function(x){return(replace(sign(x),sign(x)==0, 1))}
  # for another time: make the function robust to negatives... requires changes of < and >
  #   result = paste(ifelse(abs(sort(myvector)[2])<10^-dec, paste("<", sign0(sort(myvector)[2])*10^-dec, sep=""), format(round(sort(myvector)[2], dec), nsmall=dec)), " (", 

  result = paste(ifelse(sort(myvector)[2]<10^-dec, paste("<", 10^-dec, sep=""), format(round(sort(myvector)[2], dec), nsmall=dec)), " (", 
                 ifelse(min(myvector)<10^-dec, paste("<", 10^-dec, sep=""), format(round(min(myvector), dec), nsmall=dec)), ", ", 
                 ifelse(max(myvector)<10^-dec, paste("<", 10^-dec, sep=""), format(round(max(myvector), dec), nsmall=dec)), ")", sep="")
  return(result)
}

# ci_string_dec = function(myvector, dec){
#   result = paste(format(round(sort(myvector)[2], dec), nsmall=dec), " (", 
#                  format(round(min(myvector), dec), nsmall=dec), ", ", 
#                  format(round(max(myvector), dec), nsmall=dec), ")", sep="")
#   return(result)
# }
