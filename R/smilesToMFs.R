#' Convert SMILES to Molecular Formulas
#'
#' @param smile a character string of a valid SMILE
#' @return a \code{tibble} of Molecular Formulas
#'
#' @export
#' @examples
#' smilesToMFs('C[C@@H](C(=O)O)N')

smilesToMFs <- function(smiles)
{
  res <- map(smiles,~{
    smipar <- rcdk::parse.smiles(.)[[1]]
    
    rcdk::convert.implicit.to.explicit(smipar)
    mf <- rcdk::get.mol2formula(smipar, charge = 0)
    
    tibble(MF = mf@string)
  })
  names(res) <- smiles
  res <- res %>%
    bind_rows(.id = 'SMILE')
  return(res)
}