sim_fun_MK_RNplus<- function(pars, errors, ..., nsim = 1500){
  
  kappa_r <- pars[3]

  
  kappas <- rgamma(nsim, shape = pars[1] / pars[2], scale = pars[2])
  
  pchoose <- vector("numeric", 181)
  kc <- vector("numeric",181)
  
  for (i in seq_len(nsim)) {
    
    kc <- sqrt(kappa_r^2 + kappas[i]^2 + 2 * kappa_r*kappas[i]*cos(base_radians))
    
    pchoose <- pchoose + 
      ((besselI(kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(kappas[i],0,expon.scaled = TRUE) * 
             besselI(kappa_r,0,expon.scaled = TRUE))) *
         exp(kc - (kappas[i] + kappa_r)))
    
  }
  
  pchoose <- (pchoose / sum(pchoose)) * 0.5
  sum(log(pchoose[errors]))
}

Rcpp::sourceCpp("Rcpp/csim_fun_MK_RNplus.cpp")
