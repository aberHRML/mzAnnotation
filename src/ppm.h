#include <Rcpp.h>
using namespace Rcpp;

//' Calculate a ppm error range
//' @description Calculate the upper and lower parts per million error boundaries for a given *m/z*.
//' @param mz the *m/z* for which to calculate the error range
//' @param ppm the parts per million
//' @return A list containing the lower and upper  error range limits.
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
