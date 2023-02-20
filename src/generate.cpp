#include "elements.h"
#include "ppm.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' Calculate suitable elemental frequency ranges
//' @description Calculate elemental frequency ranges for a given mass which are suitable for molecular formula generation.
//' @param mass molecular mass
//' @return A list of minimum and maximum frequencies for each element. 
//' @examples
//' suitableElementRanges(342.11621)
//' @export
// [[Rcpp::export]]

List suitableElementRanges (float mass) {
  
  float c = round(mass/12);
  float h = round(c * 2);
  float no = round(c / 2);
  float ps = round(c / 4);
  
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

float collect_first(NumericVector vector) {
  return(vector[0]);
}

float collect_second(NumericVector vector) {
  return(vector[1]);
}

float calc_mass(vector<int> counts ,
                NumericVector element_masses,
                int charge, 
                double electron){
  int n_elements = element_masses.length();
  int i;
  float sum = 0;
  for (i=0; i < n_elements; i++)
    sum += element_masses[i] * counts[i];
  
  float mass = sum - (charge * electron);
  
  return(mass);
}

// [[Rcpp::export]]
DataFrame generate (double measured_mass,
                    double ppm,
                    int charge,
                    List element_ranges) {

  // Prepare elemental information
  DataFrame element_info = elements();
  NumericVector element_indexes = {0,2,16,19,22,23};
    
  CharacterVector all_symbols = element_info["Element"];
  CharacterVector element_symbols = all_symbols[element_indexes];
  
  NumericVector all_masses = element_info["AtomicMass"];
  NumericVector element_masses = all_masses[element_indexes];
  
  double electron = all_masses[31];
  
  // Calculate limits
  List limits = ppmRange(measured_mass,ppm);
  float limits_lower = limits["lower"];
  float limits_upper = limits["upper"];
  
  NumericVector min = sapply(element_ranges,collect_first);
  min.names() = element_ranges.names();
    
  NumericVector max = sapply(element_ranges,collect_second);
  max.names() = element_ranges.names();
  
  // Declare loop counters
  float mass = 0;
  vector<int> counts (element_indexes.length());
  
  NumericVector save (element_indexes.length());
  save.names() = element_ranges.names();
  
  long long counter = 0;
  
  // Declare output vectors
  CharacterVector generated_formulas;
  NumericVector generated_masses;
  
  // Begin loops
  counts[5] = min["S"] - 1; save["S"] = counts[5];
  while (counts[5]++ < max["S"]) { 
    counts[4] = min["P"] - 1; save["P"] = counts[4];
    while (counts[4]++ < max["P"]) { 
      counts[3] = min["O"] - 1; save["O"] = counts[3];
      while (counts[3]++ < max["O"]) { 
        counts[2] = min["N"] - 1; save["N"] = counts[2];
        while (counts[2]++ < max["N"]) {
          counts[0] = min["C"] - 1; save["C"] = counts[2];
          while (counts[0]++ < max["C"]) {
            counts[1] = min["H"] - 1; save["J"] = counts[1];
            while (counts[1]++ < max["H"]) { 
              mass = calc_mass(counts,element_masses,charge,electron);
              
              // Break loop if upper limit reached
              if (mass > limits_upper)  break;
              
              // Check the mass is between the limits
              if ((mass >= limits_lower) && (mass <= limits_upper)){ 
                
                // Create the molecular formula string
                long i;
                String mf;
                
                for (i = 0; i < element_masses.length(); i++){			
                  if (counts[i] > 0) {
                    mf += element_symbols[i];
                    
                    if (counts[i] > 1) {
                      String count = counts[i];
                      mf += count;
                    }
                  }
                }
                generated_formulas.insert(counter,mf);
                generated_masses.insert(counter,mass);
                counter++;
              }
            }
            if ((mass >= limits_lower) && (save["H"] == counts[1]-1)) break;
          }
          if ((mass >= limits_lower) && (save["C"] == counts[0]-1)) break;
        }
        if ((mass >= limits_lower) && (save["N"] == counts[2]-1)) break;
      }
      if ((mass >= limits_lower) && (save["O"] == counts[3]-1)) break;
    }
    if ((mass >= limits_lower) && (save["P"] == counts[4]-1)) break;
  }
  
  DataFrame results = DataFrame::create(
    Named("MF") = generated_formulas,
    Named("Mass") = generated_masses
  );
  
  return(results);
}
