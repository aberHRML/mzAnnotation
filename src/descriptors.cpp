#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>
#include <openbabel/descriptor.h>
using namespace OpenBabel;

#include <Rcpp.h>
using namespace Rcpp;

//' SMARTS substructure search
//' @description SMARTS substructure searching for SMILES.
//' @param smile a valid SMILE
//' @param smart a valid SMARTS symbol
//' @examples
//' smartsSearch(amino_acids$SMILES[1],"[OX2H]")
//' @export
// [[Rcpp::export]]
int smartsSearch(std::string smile,std::string smart){

  OBMol mol;
  OBConversion conv;

  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);

  OBSmartsPattern smarts;
  std::vector<std::vector<int> > maplist;

  smarts.Init(smart);
  smarts.Match(mol);
  maplist = smarts.GetMapList();
  int match = maplist.size();
  return match;
}

// [[Rcpp::export]]
double descriptor(std::string smile,const char* desc){
  
  double res;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  OBDescriptor* descType = OBDescriptor::FindType(desc);
  if(descType)
    res = descType->Predict(&mol, NULL);

  return res;
}
