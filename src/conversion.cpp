#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
using namespace OpenBabel;

#include <Rcpp.h>
using namespace Rcpp;

//' convert
//' @description convert between SMILES and Inchi and to InchiKey
//' @param input a valid SMILE or Inchi
//' @param inputType either "smiles" or "inchi", denoting the input type
//' @param outputType either "smiles", "inchi" or "inchikey", denoting the output type
//' @examples
//' convert(aminoAcids$SMILE[1],'smiles','inchi')
//' @export
// [[Rcpp::export]]
std::string convert(std::string input,const char* inputType,const char* outputType){
  
  std::string res = "NA";
  OBMol mol;
  OBConversion conv;
  
  conv.SetInAndOutFormats(inputType,outputType);

  if(conv.ReadString(&mol, input)) {
    res = conv.WriteString(&mol, true);
  }
  
  return res;
}

//' smileToMF
//' @description convert a smile to a molecular formula
//' @param smile a valid SMILE
//' @examples
//' smileToMF(aminoAcids$SMILE[1])
//' @export
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

//' smileToAccurateMass
//' @description convert a smile to an accurate mass
//' @param smile a valid SMILE
//' @examples
//' smileToAccurateMass(aminoAcids$SMILE[1])
//' @export
// [[Rcpp::export]]
double smileToAccurateMass(std::string smile){
  double res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  res = mol.GetExactMass();
  return res;
}