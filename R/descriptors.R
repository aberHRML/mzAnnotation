#' convert SMILE to a series of molecular descriptors
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
#' @importFrom ChemmineOB forEachMol smartsSearch_OB prop_OB
#' @importFrom parallel makeCluster parLapply stopCluster
#' @export
#' @examples
#' data(aminoAcids)
#' descriptors(aminoAcids$SMILE)

descriptors <- function(smiles,nCores = 1,clusterType = 'FORK',verbose = T){
  if (verbose == T) {
    cat(length(smiles),'SMILEs\n')
  }
  
  getProp <- function(smile){
    forEachMol("SMILES",smile,identity) %>%
      prop_OB() %>%
      as_tibble()
  }
  
  if (verbose == T) {
    cat('\nGenerating properties')
  }
  
  if (nCores > 1) {
    clus <- makeCluster(nCores,clusterType)
    prop <- parLapply(clus,smiles,getProp)
    stopCluster(clus)
  } else {
    prop <-  map(smiles,getProp)  
  }
  
  prop <- prop %>%
    bind_rows()
  
  if (verbose == T) {
    cat('\nGenerating functional descriptors')
  }
  
  Fgroups <- tibble(Name = c('NHH',
                             'OH',
                             'COOH',
                             'COO'),
                    String = c("[NX3;H2]",
                               "[OX2H]",
                               "[CX3](=O)[OX2H1]",
                               "[CX3](=O)[OX1H0-]")
  )
  
  
  groups <- Fgroups %>%
    split(1:nrow(Fgroups)) %>%
    map(~{
      cat('\n',.$Name)
      string <- .$String
      g <- map_int(smiles,~{
        forEachMol("SMILES",.,identity) %>%
          smartsSearch_OB("[NX3;H2]")
      }) %>%
        as_tibble()
      names(g) <- .$Name
      return(g)
    }) %>%
    bind_cols()
  
  desc <- bind_cols(prop,groups)
  
  return(desc)
}