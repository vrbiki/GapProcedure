#include <Rcpp.h>
using namespace Rcpp;

//  This function calculates the the empirical base frequencies
//  using the average value over all sequences. In other
//  words, it finds the proportion of A, G, C, and Ts
//  observed over the entire data set
//
//' Calculates the empirical base frequencies
//'
//' @param x A character vector of the sequence data.
//' @export
// [[Rcpp::export]]
List pifreq(CharacterVector x) {
  int n = x.size();
  Rcpp::IntegerVector out(4);
  int denom = 0;

  out[0] = 0;
  out[1] = 0;
  out[2] = 0;
  out[3] = 0;
  denom = x.size();
  for(int i = 0; i < n; ++i) {
    if (x[i]=="a"){
      out[0] += 1;
    } else if (x[i]=="g"){
      out[1] += 1;
    } else if(x[i]=="c") {
      out[2] += 1;
    } else if (x[i]=="t"){
      out[3] += 1;
    } else if (x[i]=="r"){
      out[0] += 1;
      out[1] += 1;
      denom += 1;
    } else if (x[i]=="m"){
      out[0] += 1;
      out[2] += 1;
      denom += 1;
    } else if (x[i]=="w"){
      out[0] += 1;
      out[3] += 1;
      denom += 1;
    } else if (x[i]=="s"){
      out[1] += 1;
      out[2] += 1;
      denom += 1;
    } else if (x[i]=="k"){
      out[1] += 1;
      out[3] += 1;
      denom += 1;
    } else if (x[i]=="y"){
      out[2] += 1;
      out[3] += 1;
      denom += 1;
    } else if (x[i]=="v"){
      out[0] += 1;
      out[1] += 1;
      out[2] += 1;
      denom += 2;
    } else if (x[i]=="h"){
      out[0] += 1;
      out[2] += 1;
      out[3] += 1;
      denom += 2;
    } else if (x[i]=="d"){
      out[0] += 1;
      out[1] += 1;
      out[3] += 1;
      denom += 2;
    } else if (x[i]=="b"){
      out[1] += 1;
      out[2] += 1;
      out[3] += 1;
      denom += 2;
    } else {
      out[0] += 1;
      out[1] += 1;
      out[2] += 1;
      out[3] += 1;
      denom += 3;
    }
  }

  List z  = List::create(out, denom);
  return z;
}
