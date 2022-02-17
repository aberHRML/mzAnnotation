#include <Rcpp.h>
using namespace Rcpp;

//' Calculate suitable elemental frequency ranges
//' @description Calculate elemental frequency ranges for a given mass which are suitable for molecular formula generation.
//' @param mass accurate mass
//' @return named numeric vector of element frequencies
//' @examples
//' suitableElementRanges(342.11621)
//' @export
// [[Rcpp::export]]

List suitableElementRanges (double mass) {
  
  double c = round(mass/12);
  double h = round(c * 2);
  double no = round(c / 2);
  double ps = round(c / 4);
  
  NumericVector carb = {0,c};
  NumericVector Hs = {0,h};
  NumericVector NO = {0,no};
  NumericVector PS = {0,ps};
  
  List maxi = List::create(Named("C") = carb,
                           Named("H") = Hs,
                           Named("N") = NO,
                           Named("O") = NO,
                           Named("P") = PS,
                           Named("S") = PS);
  
  return(maxi);
}
