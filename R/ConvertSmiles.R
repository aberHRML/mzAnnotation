#' Convert SMILES to Molecular Formula's and Exact Masses
#'
#' @param smile a character string of a valid SMILE
#' @return a \code{tibble} of Molecular Formula and Exact Mass
#'
#' @export
#' @examples
#' ConvertSmile('C[C@@H](C(=O)O)N')

ConvertSmiles <- function(smile)
{
  smipar <- rcdk::parse.smiles(smile)[[1]]
  
  rcdk::convert.implicit.to.explicit(smipar)
  mf <- rcdk::get.mol2formula(smipar, charge = 0)
  
  mf_tib <- tibble::tibble(MF = mf@string, EXACT_MASS = mf@mass)
  
  return(mf_tib)
}