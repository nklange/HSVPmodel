#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector csim_fun_MK_RNplus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  double kappa;
  double kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians;
  
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {
      
      kappa = R::rgamma(pars[0]/pars[1], pars[1]);
      
      kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(baseradians[radian]));
      
      radians = radians + (R::bessel_i(kc,0,2) / 
        (2*PI*R::bessel_i(kappa,0,2) * 
        R::bessel_i(pars[2],0,2))) *
        exp(kc - (kappa + pars[2]));
      
    }
    
    out[radian] = radians / x;
    radians = 0;
  }
  
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
