#' Find parameters from a Gamma distribution using quantiles
find_gamma_pars = function(lci, hci){
  gam_find = function(alphas){
    return(abs(qgamma(0.025, alphas, 1)/qgamma(1-0.025, alphas, 1) - lci/hci))
  }
  final_alpha = optimize(gam_find, c(0.000001, 200000)) #minimizes
  final_beta = lci/qgamma(0.025, final_alpha$minimum, 1)
  return(c(final_alpha$minimum, final_beta))
}

find_gamma_flex = function(lci, hci, bound1, bound2){
  gam_find = function(alphas){
    return(abs(qgamma(bound1, alphas, 1)/qgamma(bound2, alphas, 1) - lci/hci))
  }
  final_alpha = optimize(gam_find, c(0.000001, 200000)) #minimizes
  final_beta = lci/qgamma(bound1, final_alpha$minimum, 1)
  return(c(final_alpha$minimum, final_beta))
}