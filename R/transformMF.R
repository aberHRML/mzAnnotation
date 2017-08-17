#' transformMF
#' @description transform a molecular formula
#' @param MF molecular formula to transform
#' @param transformation transformation to apply
#' @param transformations transformations table containing available transformations rules. Defaults to \code{Transformations}.
#' @details \code{NA} will be returned if \code{MF} cannot be transformed.
#' @examples 
#' transformMF('C4H5O5')
#' transformMF('C4H5N',transformation = 'M - [OH] + [NH2]')
#' @importFrom tidyr gather
#' @importFrom stringr str_c str_replace
#' @export

transformMF <- function(MF, transformation = 'M - [O] + NH2]', transformations = mzAnnotation::Transformations){
  if (!is.na(transformation)) {
    elements <- c('C','H','O','N','P','S')
    
    transformation <- filter(transformations,`MF Change` == transformation) %>%
      select(C:S)
    
    MF <- makeup(MF)
    
    if (length(which(!(elements %in% names(MF)))) > 0) {
      MF <- c(MF,rep(0,length(which(!(elements %in% names(MF))))))
      names(MF)[names(MF) == ''] <- elements[!(elements %in% names(MF))]
    }
    
    MF <- MF[order(names(MF))]
    MF <- MF + transformation
    MF <- gather(MF,'Element','Frequency') 
    
    if (T %in% (MF$Frequency < 0)) {
      MF <- NA
    } else {
      MF <- MF %>%
        filter(Frequency > 0)
      
      MF$Frequency[MF$Frequency == 1] <- ''
      
      MF <- str_c(MF$Element,MF$Frequency) %>%
        str_c(collapse = '') 
    }
  }
  return(MF)
}