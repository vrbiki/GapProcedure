#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix AdjJC69(IntegerMatrix nDNAmat, NumericMatrix jcmat) {
  int nrows = nDNAmat.nrow();
  int ncols = nDNAmat.ncol();
  int a, b;

  // allocate the matrix we will return
  NumericMatrix rmat(nrows, nrows);

  for (int i = 0; i < nrows; i++) {
    for (int j = i+1; j < nrows; j++) {

      double count = 0.0;

      for (int s = 0; s < ncols; s++){

        if (nDNAmat(i,s)!=nDNAmat(j,s)){

          a = nDNAmat(i,s);
          b = nDNAmat(j,s);
          count += jcmat(a,b);

        } else {

          if (nDNAmat(i,s)>3){ // if its not AA, GG, CC, or TT
            a = nDNAmat(i,s);
            b = nDNAmat(j,s);
            count += jcmat(a,b);
          }

        }


      } // end s loop

      rmat(i,j) = rmat(j,i) =  (double) -.75*log(1 - 4*count/(3*ncols));

    }// end j loop
  }// end i loop
  return rmat;
}


