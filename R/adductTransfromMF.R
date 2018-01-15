#' adductTransfromMF
#' @description adduct transform a molecular formula
#' @param MF molecular formula to transform
#' @param adduct adduct to use for transformation
#' @examples 
#' adductTransfromMF('C6H12O6','[M+H]1+')
#' @export

adductTransformMF <- function(MF,adduct){

  tMF <- function(freq,trans,expres){
    if (!is.na(trans)) {
     expAt <- trans %>% makeup()
     if (F %in% (names(expAt) %in% names(freq))) {
       tmp <- rep(0,length(which(!(names(expAt) %in% names(freq)))))
       names(tmp) <- names(expAt)[!(names(expAt) %in% names(freq))]
       freq <- c(freq,tmp)
       freq <- freq[order(names(freq))]
     }
     if (expres == '+') {
       freq[names(expAt)] <- freq[names(expAt)] + expAt
     }
     if (expres == '-') {
       freq[names(expAt)] <- freq[names(expAt)] - expAt
     }
    }
    return(freq)
  }
  
  adductRule <- Adducts %>% 
    filter(Name == adduct)
  
  freq <- MF %>% makeup()
  
  freq <- freq %>% 
    tMF(adductRule$RemAt,'-') %>%
    tMF(adductRule$AddEx,'+') %>%
    tMF(adductRule$RemEx,'-')
  
  if (T %in% (freq < 0)) {
    stop('Not possible!')
  }
  
  freq <- freq[!(freq == 0)]
  freq[freq == 1] <- ''
  freq <- str_c(names(freq),freq) %>%
    str_c(collapse = '')
  
  return(freq)
  
}
