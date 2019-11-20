int_fun_MK_RNplus <- function(x, pars, radian) {
  
  out <- vector("numeric", length(x))
  
  for (i in seq_along(x)) {
    
    kappa <- x[i]
    
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

int_fun_J_RNplus <- function(x, pars, radian) {
  out <- vector("numeric", length(x))
  
  for (i in seq_along(x)) {
    
    kappa <- KappafromJ(x[i])
    kc <- sqrt(pars[3]^2 + kappa^2 + 2 * pars[3]* kappa*cos(radian))
    
    out[i] <- dgamma(x[i], shape = pars[1]/pars[2], scale = pars[2]) * 
      ((besselI(x = kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(x = kappa,0,expon.scaled = TRUE) * 
             besselI(x = pars[3],0,expon.scaled = TRUE))) * 
         exp(kc - (kappa + pars[3])))
    
  }
  out
  
}

Rcpp::sourceCpp("Rcpp/cint_fun_J_RNplus.cpp")

int_fun_J_RNminus <- function(x, pars, radian) {
  out <- vector("numeric", length(x))
  
  for (i in seq_along(x)) {
    
    kappa <- KappafromJ(x[i])
   
    out[i] <- dgamma(x[i], shape = pars[1]/pars[2], scale = pars[2]) * 
      (1 / (2 * pi * besselI(x =kappa, nu = 0, expon.scaled = TRUE)) *
         (exp(cos(radian) - 1)) ^ kappa)
  }
  out
  
}
Rcpp::sourceCpp("Rcpp/cint_fun_J_RNminus.cpp")

int_fun_MK_RNminus <- function(x, pars, radian) {
  out <- vector("numeric", length(x))
  
  for (i in seq_along(x)) {
  
    
    out[i] <- dgamma(x[i], shape = pars[1]/pars[2], scale = pars[2]) * 
      (1 / (2 * pi * besselI(x =x[i], nu = 0, expon.scaled = TRUE)) *
         (exp(cos(radian) - 1)) ^ x[i])
  }
  out
  
}
Rcpp::sourceCpp("Rcpp/cint_fun_J_RNminus.cpp")