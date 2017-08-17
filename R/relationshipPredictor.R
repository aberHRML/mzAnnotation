#' m/z relationship prediction
#' @description adduct, isotope and biotransfromation prediction.
#' @param mz numeric \code{vector} of accurate m/z
#' @param limit limit of deviation for thresholding associations. Defaults to 0.001
#' @param adducts character vector of adducts to use. If \code{NULL} all available adducts will be used.
#' @param isotopes character vector of isotopes to use.
#' @param transformations character vector of transformations to use
#' @author Jasen Finch
#' @export
#' @importFrom utils combn
#' @importFrom dplyr left_join contains
#' @importFrom stringr str_c
#' @examples 
#' relationshipPredictor(c(132.03023,168.00691))

relationshipPredictor <- function(mz, limit = 0.001, adducts = c("[M-H]1-","[M+Cl]1-","[M+K-2H]1-"), isotopes = NULL, transformations = NULL, adductTable = mzAnnotation::Adducts, isotopeTable = mzAnnotation::Isotopes, transformationTable = mzAnnotation::Transformations){
  if (is.null(adducts)) {
    adducts <- adductTable$Name
  }
  isotopes <- c('NA',isotopes)
  transformations <- c('NA',transformations)
  
  combinations <- combn(mz,2)
  
  combinations <- apply(combinations,2,function(x){
    M1 <- calculateMs(x[1],adducts,isotopes,transformations,adductTable,isotopeTable,transformationTable)
    M2 <- calculateMs(x[2],adducts,isotopes,transformations,adductTable,isotopeTable,transformationTable)
    
    coms <- expand.grid(M1$ID,M2$ID)
    colnames(coms) <- c('ID1','ID2')
    coms <- left_join(coms,M1,by = c('ID1' = 'ID'))
    colnames(coms)[3:ncol(coms)] <- str_c(colnames(coms)[3:ncol(coms)],'1')
    coms <- left_join(coms,M2,by = c('ID2' = 'ID'))
    colnames(coms)[8:ncol(coms)] <- str_c(colnames(coms)[8:ncol(coms)],'2')
    
    coms <- mutate(coms,Error = abs(M1 - M2)) %>%
      filter(Error <= limit) %>%
      select(contains('m/z'),contains('Adduct'),contains('Isotope'),contains('Transformation'),Error)
    return(coms)
  })
  combinations <- bind_rows(combinations)
  return(combinations)
}