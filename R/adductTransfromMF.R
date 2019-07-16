#' adductTransformMF
#' @description adduct transform a molecular formula.
#' @param MF molecular formula to transform.
#' @param adduct adduct to use for transformation.
#' @param Adducts Adduct information table to use. Defaults to \code{mzAnnotation::Adducts}.
#' @examples 
#' adductTransformMF('C6H12O6','[M+H]1+')
#' @export

adductTransformMF <- function(MF,adduct,Adducts = adducts()){

  tMF <- function(freq,trans,expres){
    if (!is.na(trans)) {
     expAt <- trans %>% count.elements()
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
  
  freq <- MF %>% count.elements()
  
  freq <- freq %>% 
    {. * adductRule$xM} %>%
    tMF(adductRule$RemAt,'-') %>%
    tMF(adductRule$AddEx,'+') %>%
    tMF(adductRule$RemEx,'-')
  
  if (T %in% (freq < 0)) {
    freq <- ''
  } else {
    ch <- freq[names(freq) %in% c('C','H')]
    freq <- freq[!(names(freq) %in% c('C','H'))]
    freq <- c(ch,freq)
    freq <- freq[!(freq == 0)]
    freq[freq == 1] <- ''
    freq <- str_c(names(freq),freq) %>%
      str_c(collapse = '') 
  }
  
  return(freq)
  
}
