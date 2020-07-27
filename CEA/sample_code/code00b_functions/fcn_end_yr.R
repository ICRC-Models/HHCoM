end_yr = function(df, yrs_active_off, thresh_end){
  # thresh_end is the number of cases that would warrant shutting active surveillance down
  # yrs_active_off is the number of years in the sequence before we stop 
  # df sequence of years with the counts, rows must be labeled with calendar year
  
  # how to ask for year of declaration? 
  # with rle (run length encoding), figure out the sequences
  # when is there a value of 0 that occurs more than N times
  
  # which(tmp_rle$values==0 & tmp_rle$length>3)
  # ask for the year before the sequence of zeros (so subtract 1 to value above).
  # names(tmp_rle$lengths)[which(tmp_rle$values==0 & tmp_rle$length>3)-1]
  # then add N+1 for the year when AS is seized
  # yr_declaration = as.numeric(names(tmp_rle$lengths)[which(tmp_rle$values==0 & tmp_rle$length>3)-1])+4
  
  
  tmp_rle = rle(df)
  tmp_rle_ofinterest = tmp_rle$values==thresh_end & tmp_rle$length>yrs_active_off
  tmp_abort = ifelse(length(which(tmp_rle_ofinterest))>0, 
                     which(tmp_rle_ofinterest), 
                     length(tmp_rle$values))
  
  if(length(tmp_rle$values)>1){
  return(min(as.numeric(names(tmp_rle$lengths)[tmp_abort-1])+(yrs_active_off), max(as.numeric(names(df)))))
  } else if(tmp_rle$values==0 & length(tmp_rle$values)==1) {
  return(as.numeric(names(df)[yrs_active_off+1])) 
  } else {
  return(max(as.numeric(names(df))))
  }
}


