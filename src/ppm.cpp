#include <Rcpp.h>
using namespace Rcpp;

//' Calculate a PPM error range
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

//' Calculate PPM error
//' @description Calculate ppmError between a measured and theoretical m/z
//' @param measured measured m/z
//' @param theoretical theoretical m/z
//' @examples ppmError(118.08626,118.08647)
//' @export
// [[Rcpp::export]]

double ppmError(double measured,double theoretical) {
  double error = (measured - theoretical) / theoretical * pow(10,6);
  return(error);
}
