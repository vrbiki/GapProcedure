#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix isNeighbours(IntegerMatrix ordermat, IntegerVector intCutoff) {
  int N = ordermat.nrow();
  IntegerMatrix rmat(N, N);
  int x = 0, y = 0;


  for(int k = 0; k < N; k++) {
    int cutoff = intCutoff[k];

    // Rcout << "k =  " << k << " Cutoff " << cutoff << "****************** "<< std::endl;

    for (int i = 0; i < cutoff; i++){
      x = ordermat(k,i);

      for(int j = 0; j < cutoff; j++) {

        y = ordermat(k,j);

        // Rcout << " x=" << x << " y=" << y <<std::endl;

        if (x==(k+1)){
          // Rcout << " rmat(" << k+1 << "," << y << ") gets 1" <<std::endl;
          rmat(k, y-1) = 1;
        }
        if (y==(k+1)) {
          // Rcout << " rmat(" << k+1 << "," << x << ") gets 1" <<std::endl;
          rmat(k, x-1) = 1;
        }

      } // end j loop

    } // end k loop

  } // end i loop

  return(rmat);

}
