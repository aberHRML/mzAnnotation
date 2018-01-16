#' Extracting count data for functional groups
#'
#' @param x a character string of a valid SMILE
#' @return a \code{tibble} of functional group counts
#' 
#' @export
#' @examples 
#' functionalGroups('CN1C=C(N=C1)CC(C(=O)O)N')
#' 

functionalGroups <- function(x)
  {
  
  molRefs = ChemmineOB::forEachMol("SMILES",x,identity)
  
  NH2 <- ChemmineOB::smartsSearch_OB(molRefs,"[NX3;H2]")
  OH <- ChemmineOB::smartsSearch_OB(molRefs,"[OX2H]")
  COOH <- ChemmineOB::smartsSearch_OB(molRefs,"[CX3](=O)[OX2H1]")
  COO <- ChemmineOB::smartsSearch_OB(molRefs,"[CX3](=O)[OX1H0-]")
  
  return(tibble::tibble(NH2, OH, COOH, COO))
  }


