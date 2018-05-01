#' Convert SMILES to Accurate Masses
#'
#' @param smile a character string of a valid SMILE
#' @return a \code{tibble} of Accurate Masses
#'
#' @export
#' @examples
#' smilesToAccurateMasses('C[C@@H](C(=O)O)N')

smilesToAccurateMasses <- function(smiles)
{
  res <- map(smiles,~{
    smipar <- rcdk::parse.smiles(.)[[1]]
    
    rcdk::convert.implicit.to.explicit(smipar)
    mf <- rcdk::get.mol2formula(smipar, charge = 0)
    
    tibble(AccurateMass = mf@mass)
  })
  names(res) <- smiles
  res <- res %>%
    bind_rows(.id = 'SMILE')
 return(res)
}