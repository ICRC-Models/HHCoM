## A function to turn an RGB code into the hex codes used by most of R.


colorvector = function(x){
  # x can be given as a vector of three elements, which yields one color.
  # In that case, then the function is equivalent to the colorspace::RGB function
  # the bigger reason to create this function was to be able to feed rgb a dataframe or matrix of colors.
  if(is.null(dim(x))){
    out = rgb(x["R"],x["G"], x["B"], max=255)
  } else {
    out = c()
    for(i in 1:dim(x)[1]){out[i] = rgb(x[i,"R"],x[i,"G"], x[i,"B"], max=255)}
  }
  return(out)
}

# test: colorvector(c("R"=1,"G"=2,"B"=3))
# test: colorvector(data.frame(R=1:3, G=4:6, B=7:9))
