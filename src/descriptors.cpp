#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>
using namespace OpenBabel;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int smartsSearch(std::string smile,std::string smart){
  
  int match;
  OBMol mol;
  OBConversion conv;
  
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smile);
  
  OBSmartsPattern smarts;
  std::vector<std::vector<int> > maplist;
  
  smarts.Init(smart);
  smarts.Match(mol);
  maplist = smarts.GetMapList();
  match = maplist.size();
  return match;
}

