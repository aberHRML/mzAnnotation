#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
using namespace OpenBabel;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string cnvrt(std::string input,const char* inputType,const char* outputType){
  
  std::string res = "NA";
  OBMol mol;
  OBConversion conv;
  
  conv.SetInAndOutFormats(inputType,outputType);

  if(conv.ReadString(&mol, input)) {
    res = conv.WriteString(&mol, true);
  }
  
  return res;
}

//' Convert SMILES to molecular formula
//' @description convert a smile to a molecular formula
//' @param smile a valid SMILE
//' @examples
//' smilesToMF(amino_acids$SMILES[1])
//' @export
// [[Rcpp::export]]
std::string smilesToMF(std::string smile){
  std::string res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  res = mol.GetFormula();
  return res;
}

//' Convert SMILES to accurate mass
//' @description convert a smile to an accurate mass
//' @param smile a valid SMILE
//' @examples
//' smilesToAccurateMass(amino_acids$SMILES[1])
//' @export
// [[Rcpp::export]]
double smilesToAccurateMass(std::string smile){
  double res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  res = mol.GetExactMass();
  return res;
}
