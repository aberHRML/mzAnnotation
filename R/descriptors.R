#' Convert SMILE to a series of molecular descriptors
#'
#' @param smile a valid SMILE
#' @return a \code{tibble} of the following;
#' \itemize{
#'     \item{HBD} The number of Hydrogen Bond Donors
#'     \item{HBA} The number of Hydrogen Bond Acceptors
#'     \item{AcidGroups} The number of Acidic Groups
#'     \item{BaseGroups} The number of Basic Groups
#'     \item{TotalCharge} The total overall charge
#' }
#'
#' @export
#' @examples
#' descriptors('CN1C=C(N=C1)CC(C(=O)O)N')

descriptors <- function(smile)
{
  smi <- rcdk::parse.smiles(smile)[[1]]
  
  rcdk::convert.implicit.to.explicit(smi)
  
  HBD <-
    rcdk::eval.desc(
      smi,
      'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor'
    )[[1]]
  HBA <-
    rcdk::eval.desc(
      smi,
      'org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor'
    )[[1]]
  
  AcidGroups <-
    rcdk::eval.desc(
      smi,
      'org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor'
    )[[1]]
  
  BaseGroups <-
    rcdk::eval.desc(
      smi,
      'org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor'
    )[[1]]
  
  TotalCharge <- rcdk::get.total.charge(smi)
  
  desc_tib <- tibble::tibble(HBD, HBA, AcidGroups, BaseGroups, TotalCharge)
  
  return(desc_tib)
}