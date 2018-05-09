#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
using namespace OpenBabel;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string convert(std::string input,const char* inputType,const char* outputType){
  
  std::string res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInAndOutFormats(inputType,outputType);
  conv.ReadString(&mol, input);
  
  res = conv.WriteString(&mol, true);
  
  return res;
}

// [[Rcpp::export]]
std::string smileToMF(std::string smile){
  std::string res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  res = mol.GetFormula();
  return res;
}
