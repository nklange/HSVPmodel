library("tidyverse")
library("GA")
library("circular")
library("parallel")
library("microbenchmark")


base_radians <- circular((0:180) * 2 * pi/360)

kappa_map <- c(seq(0, 10, length.out = 250),seq(10.001,1e4,length.out = 250))
J_map <- kappa_map*besselI(kappa_map,1,expon.scaled = TRUE)/besselI(kappa_map,0,expon.scaled = TRUE)
mapJkappa <- approxfun(J_map,kappa_map, yleft = 0)

KappafromJ <- function(J) {
  out <- mapJkappa(J)
  ifelse(is.na(out), J, out)
}


# simJRNplus --------------------
sim_j_plus_r <- function(x, pars, baseradians) {
  
  kappa_r <- pars[3]
  rJ <- rgamma(x, shape = pars[1]/pars[2], scale = pars[2])
  kappas <- KappafromJ(rJ)
  
  pchoose <- vector("numeric", 181)
  kc <- vector("numeric",181)
  
  for (i in seq_len(x)) {
    
    kc <- sqrt(kappa_r^2 + kappas[i]^2 + 2 * kappa_r*kappas[i]*cos(baseradians))
    
    pchoose <- pchoose + 
      ((besselI(kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(kappas[i],0,expon.scaled = TRUE) * 
             besselI(kappa_r,0,expon.scaled = TRUE))) *
         exp(kc - (kappas[i] + kappa_r)))
    
  }
  pchoose/sum(pchoose)
  
}
Rcpp::sourceCpp("Rcpp/csim_fun_J_RNplus.cpp")


microbenchmark::microbenchmark(
  sim_j_plus_r(1000, pars = c(20, 2, 3),baseradians = base_radians),
  csim_fun_J_RNplus(1000, pars = c(20, 2, 3),baseradians = base_radians)
)

#rcpp version needs normalizing of results and default to J at >kappa_map
#currently R 1: 2 Rcpp meh

# numintJRNplus ----------------------------------

int_j_plus_r <- function(x, pars, radian) {
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

microbenchmark::microbenchmark(
int_j_plus_r(c(0,1,10,100), pars = c(20, 2, 30),radian=0.5),
cint_fun_J_RNplus(c(0,1,10,100), pars = c(20, 2, 30),radian=0.5)
)

# currently R 4 : 1 Rcpp
# rcpp version needs normalizing of results and default to J at >kappa_map
# looks like numint benefits from rcpp but simulation doesn't