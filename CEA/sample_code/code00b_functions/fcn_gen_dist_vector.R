gen_dist_vector = function(desc, iter){
  # desc must constain:
  # rcode_dist # the kind of distribution desired: rnorm, rbeta, rgamma, etc
  # rcode_par1 # first parameter
  # rcode_par2 # second parameter
  if (desc$rcode_dist != "rgamma"){
    eval(parse(text=paste("dist_vector = ", desc$rcode_dist, "(", iter, "," , desc$rcode_par1, ",", desc$rcode_par2, ")", sep="")))
  } else {
    eval(parse(text=paste("dist_vector = ", desc$rcode_dist, "(", iter, "," , desc$rcode_par1, ", scale=", desc$rcode_par2, ")", sep="")))
  }
  return(dist_vector)
}