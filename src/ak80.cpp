#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix AdjK80(IntegerMatrix nDNAmat, NumericMatrix tsmat, NumericMatrix tvmat) {
  int nrows = nDNAmat.nrow();
  int ncols = nDNAmat.ncol();
  int a, b;

  // allocate the matrix we will return
  NumericMatrix rmat(nrows, nrows);

  for (int i = 0; i < nrows; i++) {
    for (int j = i+1; j < nrows; j++) {

      double countts = 0.0;
      double counttv = 0.0;

      for (int s = 0; s < ncols; s++) {

        if (nDNAmat(i,s)!=nDNAmat(j,s)){

          a = nDNAmat(i,s);
          b = nDNAmat(j,s);
          countts += tsmat(a,b);
          counttv += tvmat(a,b);

        } else {

          if (nDNAmat(i,s)>3){ // if its not AA, GG, CC, or TT
            a = nDNAmat(i,s);
            b = nDNAmat(j,s);
            countts += tsmat(a,b);
            counttv += tvmat(a,b);
          }

        }

      }// end s loop

      double term1 = -0.5*log(1-2*(countts/ncols)-(counttv/ncols));
      double term2 = .25*log(1-2*(counttv/ncols));
      rmat(i,j) = rmat(j,i) =  term1-term2;

    }// end j loop
  }// end i loop

  return rmat;
}


