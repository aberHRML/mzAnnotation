#include <Rcpp.h>
using namespace Rcpp;

//' ppmRange
//' @description Calculate the upper and lower ppm boundaries for a given m/z
//' @param mz the m/z for which to calculate the range
//' @param ppm the ppm 
//' @examples ppmRange(118.08626,5)
//' @export
// [[Rcpp::export]]

List ppmRange(double mz,double ppm) {
  double amu = mz * ppm * pow(10,-6);
  NumericVector lower = NumericVector::create(mz - amu);
  lower = round(lower,5);
  NumericVector upper = NumericVector::create(mz + amu);
  upper = round(upper,5);
  List res = List::create(lower(0),upper(0));
  res.attr("names") = CharacterVector::create("lower","upper");
  return(res);
}