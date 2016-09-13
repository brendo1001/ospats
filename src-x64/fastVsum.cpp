#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastVsum(NumericMatrix Ar) {
  int n = Ar.nrow();
  NumericMatrix mat(n, n);
  
  for(int i = 0; i < n; ++i) {
    mat(i,_) = Ar(i,0) + Ar(_,0);
    }
  return(mat);
}