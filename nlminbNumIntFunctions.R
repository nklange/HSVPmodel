int_fun_MK_RNplus <- function(x, pars, radian) {
  
  out <- vector("numeric", length(x))
  
  for (i in seq_along(x)) {
    
    kappa <- x[i]
    #kappa <- min(x[i],max(kappa_map))
    
    kc <- sqrt(pars[3]^2 + kappa^2 + 2 * pars[3]*kappa*cos(radian))
    
    out[i] <- dgamma(kappa, shape = pars[1]/pars[2], scale = pars[2]) 
    if (out[i] != 0) {
      out[i] <- out[i] * ((besselI(kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(kappa,0,expon.scaled = TRUE) * 
             besselI(pars[3],0,expon.scaled = TRUE))) *
         exp(kc - (kappa + pars[3])))
    }
    
  }
  out
}


Rcpp::sourceCpp("Rcpp/cint_fun_MK_RNplus.cpp")
