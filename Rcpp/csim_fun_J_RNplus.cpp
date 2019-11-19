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

NumericVector LinearSpacedArray(double a, double b, double N)
{
  double h = (b - a) / (N-1);
  NumericVector xs(N);
  int x;
  for (x = 0; x < N; x++) {
    
    if (x == 0) {
      xs[x] = a;
    } else {
      xs[x] = xs[x-1] + h;
    }
  }
  return xs;
  
}


NumericVector allK(NumericVector kSmall, NumericVector kLarge){
  
  NumericVector allK(500);
  int x;
  
  for (x = 0; x < 500; x++){
    
    if (x < 250){
      allK[x] = kSmall[x];
    } else {
      allK[x] = kLarge[x-250];
    }
    
  }
  
  return(allK);
}


NumericVector kSmall = LinearSpacedArray(0.001,10,250);
NumericVector kLarge = LinearSpacedArray(10.001,1e4,250);
NumericVector Kmap = allK(kSmall,kLarge);


NumericVector allJ(NumericVector Kmap){
  
  NumericVector allJ(500);
  int x;
  
  for (x = 0; x < 500; x++){
    
    allJ[x] = Kmap[x] * R::bessel_i(Kmap[x],1,2)/R::bessel_i(Kmap[x],0,2);
    
  }
  
  return(allJ);
}

NumericVector Jmap = allJ(Kmap);


double cKappaFromJ(double J){
  
  double v = J;
  NumericVector y = Kmap;
  NumericVector x = Jmap;
  
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int n = Jmap.size();
  
  
  int i, j, ij;
  
  //if(!n) return R_NaN;
  
  i = 0;
  j = n - 1;
  
  /* handle out-of-domain points */
  //if(v < x[i]) return Meth->ylow;
  //if(v > x[j]) return Meth->yhigh;
  
  /* find the correct interval by bisection */
  while(i < j - 1) { /* x[i] <= v <= x[j] */
  ij = (i + j)/2; /* i+1 <= ij <= j-1 */
  if(v < x[ij]) j = ij; else i = ij;
  /* still i < j */
  }
  /* provably have i == j-1 */
  
  /* interpolation */
  
  if(v == x[j]) return y[j];
  if(v == x[i]) return y[i];
  /* impossible: if(x[j] == x[i]) return y[i]; */
  
  //if(Meth->kind == 1) /* linear */
  return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
  // else /* 2 : constant */
  //return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0)
  //  + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
}



// [[Rcpp::export]]
NumericVector csim_fun_J_RNplus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  double kappa;
  NumericVector J;
  double kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  J = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {
      
      kappa = cKappaFromJ(J[i]);
      
      kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(baseradians[radian]));
      
      // Rcout << kc << std::endl;
      
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

// // [[Rcpp::export]]
// void rcpp_rcout(double x){
// 
//   // printing value of vector
//   Rcout << "The value of v : " << Jmap << "\n";
// }
