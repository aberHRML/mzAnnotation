#' convert SMILE to a series of molecular descriptors for PIPs
#'
#' @param smile vector of valid SMILEs
#' @return a \code{tibble} of the following;
#' \itemize{
#'     \item{HBD} The number of Hydrogen Bond Donors
#'     \item{HBA} The number of Hydrogen Bond Acceptors
#'     \item{AcidGroups} The number of Acidic Groups
#'     \item{BaseGroups} The number of Basic Groups
#'     \item{TotalCharge} The total overall charge
#' }
#'
#' @importFrom rcdk parse.smiles convert.implicit.to.explicit eval.desc get.total.charge
#' @export
#' @examples
#' generatePIPrules('CN1C=C(N=C1)CC(C(=O)O)N')

generatePIPrules <- function(smiles)
{
  
  smi <- map(smiles,~{
    parse.smiles(.) %>%
      .[[1]]
    })
  
  desc <- tibble(
    CHNOPS = '',
    HBD = map_int(smi,~{
      s <- .
      convert.implicit.to.explicit(s)
      eval.desc(
        s,
        'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor'
      )[[1]]
    }),
    HBA = map_int(smi,~{
      s <- .
      convert.implicit.to.explicit(s)
      eval.desc(
        s,
        'org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor'
      )[[1]]
    }),
    AcidGroups = map_int(smi,~{
      s <- .
      convert.implicit.to.explicit(s)
      eval.desc(
        s,
        'org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor'
      )[[1]]
    }),
    BaseGroups = map_int(smi,~{
      s <- .
      convert.implicit.to.explicit(s)
      eval.desc(
        s,
        'org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor'
      )[[1]]
    }),
    TotalCharge = map_int(smi,~{
      s <- .
      convert.implicit.to.explicit(s)
      get.total.charge(s)
      })
  ) %>%
    bind_cols(smilesToMFs(smiles) %>% select(MF),
              smilesToAccurateMasses(smiles) %>% select(AccurateMass)) %>%
    mutate(SMILE = smiles) %>%
    select(SMILE,MF,AccurateMass,CHNOPS:TotalCharge)
  
  return(desc)
}